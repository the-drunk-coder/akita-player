#include <iostream>
#include <algorithm>
#include <iterator>
#include <sndfile.hh>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include "getch.h"
#include "RtAudio.h"
#include <mutex>              
#include <condition_variable>
#include <string>
#include <functional>
#include <climits>
#include "akita_filters.h"
#include <stk/FreeVerb.h>

namespace po = boost::program_options;
namespace lfree = boost::lockfree;


// format enum
enum RWTYPES { UCHAR, SHORT, FLOAT, DOUBLE };

// custom stream to extract enum from command line
std::istream &operator>>(std::istream &in,  RWTYPES &rwtype) {
  std::string token;
  in >> token;

  boost::to_upper(token);

  if (token == "UCHAR") {
    rwtype = UCHAR;
  } else if (token == "SHORT") {
    rwtype = SHORT;
  } else if (token == "FLOAT") {
    rwtype = FLOAT;
  } else if (token == "DOUBLE") {
    rwtype = DOUBLE;
  }

  return in;
}

// map enable nicer type output ...
std::map <RWTYPES, std::string> rwtypes_strings {
  { UCHAR, "uchar (8 Bit)" },
  { SHORT, "short (16 Bit)" },
  { FLOAT, "float (32 Bit)" },
  { DOUBLE, "double (64 Bit)" }
};

// to faciliate initialization ...
std::map <RWTYPES, RtAudioFormat> rwtypes_rtaudio {
  { UCHAR, RTAUDIO_SINT8 },
  { SHORT, RTAUDIO_SINT16 },
  { FLOAT, RTAUDIO_FLOAT32 },
  { DOUBLE, RTAUDIO_FLOAT64 }
};

// the current state of the noise source
enum PSTATE { PLAY, LOOP_REC, LOOP, LOOP_SILENCE, LOOP_BLOCK, SILENCE, BLOCK };

// custom stream to extract enum from command line
std::istream &operator>>(std::istream &in, PSTATE &pstate) {
  std::string token;
  in >> token;

  boost::to_upper(token);

  // only need those two
  if (token == "PLAY") {
    pstate = PLAY;
  } else if (token == "SILENCE") {
    pstate = SILENCE;
  } 

  return in;
}

// RAW: no filter client, source goes straight to output
// PROXY: signal is re-routed through filter client
// MILD: filter client weeds out non-playable parts (NaN etc.)
enum PMODE { RAW, PROXY, MILD };

// custom stream to extract enum from command line
std::istream &operator>>(std::istream &in, PMODE &pmode) {
  std::string token;
  in >> token;

  boost::to_upper(token);

  // only need those two
  if (token == "RAW") {
    pmode = RAW;
  } else if (token == "PROXY") {
    pmode = PROXY;
  } else if (token == "MILD") {
    pmode = MILD;
  } 

  return in;
}

// commands to control audio threads 
namespace COMMAND {
  enum COMMAND { MODE_CHANGE, STATE_CHANGE, SAMPLERATE_CHANGE, FUZZINESS_CHANGE, GAIN_CHANGE, MEAN_FILTER_ON, MEAN_FILTER_OFF, FILTER_ON, FILTER_OFF, LOOP_INIT, LOOP_FINISH, LOOP_RELEASE };
}

/*
 * container for command line options
 */
struct options_container {

  // the file we're working on 
  std::string filename;

  unsigned int samplerate = 44100;
  unsigned int buffer_frames = 1024;
  
  PSTATE initial_state;
  PMODE initial_mode;
  float initial_gain;
  
  // format glitch parameters
  RWTYPES read_type;
  RWTYPES write_type;
  RWTYPES stream_type;

  // buffer glitch parameters
  float buffer_cut;
  float offset_cut;
  float sample_repeat;
  float flip_prob;
  // kill samples to make it all fuzzy !
  float fuzziness;
  
  int channel_offset;

  int mean_filter_points;

  float start;
  float end;

  bool mono;
};

// file container to load and store samples in buffer
template <typename READ_TYPE>
struct file_container {
  enum STATE{ INIT, READY, FAILED };

  STATE state = INIT;

  const int CHUNKSIZE = 1024;
  
  // fields
  READ_TYPE *file_buffer;
  std::string name;
  long int samples;
  long int frames;
  long int start_sample;
  long int end_sample;
  short channels;
  int samplerate;

  // methods
  file_container(std::string fname, float start, float end, bool mono) {
    name = fname;
    // libsndfile soundfilehandle ... 
    SndfileHandle file = SndfileHandle(name.c_str());

    channels = mono ? 1 : file.channels();
    // transfer some info
    samplerate = file.samplerate();

    frames = file.frames();

    // can't work on an empty file !
    if(frames <= 0){
      state = FAILED;
      return;
    }
    
    samples = frames * channels;

    start_sample = (float) frames * start * channels;
    end_sample = (float) frames * end * channels;
  
    file_buffer = new READ_TYPE[samples];

    int frame_chunksize = CHUNKSIZE * file.channels();

    READ_TYPE *chunk_buffer = new READ_TYPE[frame_chunksize];

    // Read file into buffer -- could this be made faster ?
    long int read_offset = 0;
    while (file.readf(chunk_buffer, CHUNKSIZE) == CHUNKSIZE) {      
      if(mono) {
	int chunk_index = 0;
	for (int i = 0; i < CHUNKSIZE; i++) {
	  file_buffer[read_offset + i] = 0;
	  for(int j = 0; j < file.channels(); j++) {
	    file_buffer[read_offset + i] += chunk_buffer[chunk_index + j] / file.channels();
	    
	  }
	  chunk_index += file.channels();
	}
	read_offset += CHUNKSIZE;
      } else {
	for (int i = 0; i < frame_chunksize; i++) {	
	  file_buffer[read_offset + i] = chunk_buffer[i];
	}
	read_offset += frame_chunksize;
      }
    }
    // not needed any longer !
    delete[] chunk_buffer;
    
    state = READY;
  }
  
  ~file_container() {

    if(state == READY){
      delete[] file_buffer;
    }
  }

  void print_file_info() {
    // display some file info ...
    std::cout << "Input file: " << name << ", " << channels << "ch, "
	      << samplerate << "Hz, " << frames << " frames.\n"
	      << std::endl;
  }

};

/*
 * Command Containers for source- and filter threads.
 */
struct source_command_container {
  COMMAND::COMMAND cmd;
  
  // parameters that could be subject to change
  PSTATE new_state;
  float new_fuzziness;
  int new_sample_repeat;
};

struct filter_command_container {
  COMMAND::COMMAND cmd;
  
  // parameters that could be subject to change
  PMODE new_mode;
  float gain;  
};

/*
 * Parameter structs to pass to noise- and filter threads.
 */

// params to pass to noise source
template <typename READ_TYPE>
struct source_params {

  // the command queue
  lfree::spsc_queue<source_command_container>* cmd_queue;

  file_container<READ_TYPE>* fc;
  
  // state ... playing, silent, blocking, looping etc
  PSTATE state = PLAY;

  // loop params
  long int loop_start = 0;
  long int loop_end = 0;

  // the current position within the file buffer
  long int offset = 0;

  // buffer glitch parameters
  float buffer_cut;
  float offset_cut;
  float sample_repeat;

  float fuzziness;
  
  source_params(options_container& opts) {
    state = opts.initial_state;
    sample_repeat = opts.sample_repeat;
    fuzziness = opts.fuzziness;
    
    fc = new file_container<READ_TYPE>(opts.filename, opts.start, opts.end, opts.mono);

    // set starting point
    offset = fc->start_sample;
    
    // default handling ...
    if(opts.offset_cut == 2){      
      offset_cut = fc->channels;
      opts.offset_cut = offset_cut;
    } else {
      offset_cut = opts.offset_cut;
    }

    if(opts.buffer_cut == 2){      
      buffer_cut = fc->channels;
      opts.buffer_cut = buffer_cut;
    } else {
      buffer_cut = opts.buffer_cut;
    }
           
    cmd_queue = new lfree::spsc_queue<source_command_container>(10);
  }

  ~source_params(){
    //std::cout << "cleaning source params" << std::endl;
    delete fc;
    delete cmd_queue;
  }
};

// params for the filter 
struct filter_params {
  lfree::spsc_queue<filter_command_container>* cmd_queue;
  filterbank* fbank; 
  mean_filterbank* m_fbank;
  stk::FreeVerb rev;
  
  PMODE mode;
  int channels;
  float gain;
  bool filter = false;
  bool mean_filter = false;

  float flippiness;
  
  filter_params (options_container& opts, int channels) {
    mode = opts.initial_mode;
    gain = opts.initial_gain;
    flippiness = opts.flip_prob;
    this->channels = channels;
    
    cmd_queue = new lfree::spsc_queue<filter_command_container>(10);
    fbank = new filterbank(channels, opts.samplerate, 100, 8000, 10);
    //rev = new fv3::nrev_f;
    //rev->setrt60(1);
    
    
    if(opts.mean_filter_points > 0){
      m_fbank = new mean_filterbank(channels, opts.mean_filter_points);
      mean_filter = true;
    } else {
      m_fbank = new mean_filterbank(channels, 13);
    }
  }

  ~filter_params () {
    //std::cout << "cleaning filter params" << std::endl;
    delete cmd_queue;
    delete fbank;
    delete m_fbank;
    //delete rev;    
  }
};

// mutex and cv to block audio thread 
std::mutex mtx;
std::condition_variable cv;

/*
 * RtAudio callback functions for source- and filter clients.
 */

// The parameterized generator callback function ...
template <typename READ_TYPE, typename WRITE_TYPE>
int source_callback(void *outputBuffer, void *inputBuffer,
                 unsigned int nBufferFrames, double streamTime,
                 RtAudioStreamStatus status, void *userData) {
  
  // get the parameter container from the user data ...
  source_params<READ_TYPE> *spar = reinterpret_cast<source_params<READ_TYPE> *>(userData);
  lfree::spsc_queue<source_command_container>* cmd_queue = spar->cmd_queue;
  file_container<READ_TYPE>* fc = spar->fc;

  WRITE_TYPE *out_buf = (WRITE_TYPE *) outputBuffer;

  // handle commands
  source_command_container scont;
  while(cmd_queue->pop(&scont)){
    if(scont.cmd == COMMAND::STATE_CHANGE){
      spar->state = scont.new_state;
    } else if (scont.cmd == COMMAND::LOOP_INIT) {
      spar->loop_start = spar->offset;
      spar->state = LOOP_REC;
    } else if (scont.cmd == COMMAND::LOOP_FINISH) {
      spar->loop_end = spar->offset;
      spar->state = LOOP;
    } else if (scont.cmd == COMMAND::FUZZINESS_CHANGE){
      spar->fuzziness = scont.new_fuzziness;
    } else if (scont.cmd == COMMAND::SAMPLERATE_CHANGE){
      spar->sample_repeat = scont.new_sample_repeat;
    }
  }

  // block source thread ... might have an interesting effect ...
  if(spar->state == BLOCK || spar->state == LOOP_BLOCK){
    std::unique_lock<std::mutex> lck(mtx);
    cv.wait(lck);
    return 0;
  }

  if (status) { std::cout << "Stream underflow detected!" << std::endl; }

  // SILENCE ! I KILL YOU !!
  if (spar->state == SILENCE || spar->state == LOOP_SILENCE) {
    for (int i = 0; i < nBufferFrames * spar->buffer_cut; i++) {
      *out_buf++ = 0;
    }
    // in this case, offset won't be modified
    return 0;
  }

  // transfer samples from file buffer to output !
  for (uint32_t i = 0; i < nBufferFrames * spar->fc->channels; i++) {   
    // loop in case we hit the file's end
    if (spar->offset + i > fc->end_sample) {
      spar->offset = fc->start_sample;
    }
    // copy samples        
    out_buf[i] = fc->file_buffer[spar->offset + i];    
  }
  
  // raw sample repetition vs channel-wise sample repetition ?
  // sample repetition
  for (uint32_t i = 0; i < (nBufferFrames * spar->fc->channels) / spar->sample_repeat; i++) {
    for(int j = 0; j < spar->sample_repeat; j++){
      out_buf[i+j] = out_buf[i];
    }
  }
  
  // buffercutting
  for (uint32_t i = nBufferFrames * spar->buffer_cut; i < nBufferFrames * spar->fc->channels; i++) {           
    out_buf[i] = 0.0;    
  }
  
  // random sample shootout
  if(spar->fuzziness > 0.0){
    for (uint32_t i = 0; i < nBufferFrames * spar->fc->channels; i++) {   
      if(rand() / (float) INT_MAX < spar->fuzziness){
	out_buf[i] = 0.0;    
      }
    }
  }
  
  // increament read offset
  spar->offset += ((float) nBufferFrames * spar->offset_cut) / spar->sample_repeat;
  if((spar->state == LOOP) && (spar->offset >= spar->loop_end)){
    spar->offset = spar->loop_start;
  }
   
  return 0;
}


// the callback function for the filter thread
int filter_callback(void *outputBuffer, void *inputBuffer,
                 unsigned int nBufferFrames, double streamTime,
	   RtAudioStreamStatus status, void *userData) {

  filter_params *fpar = reinterpret_cast<filter_params*>(userData);
  lfree::spsc_queue<filter_command_container>* cmd_queue = fpar->cmd_queue;

  // handle commands
  filter_command_container fcont;
  while(cmd_queue->pop(&fcont)){
    if(fcont.cmd == COMMAND::GAIN_CHANGE) {
      fpar->gain = fcont.gain;
    } else if(fcont.cmd == COMMAND::FILTER_ON) {
      fpar->filter = true;
    } else if(fcont.cmd == COMMAND::FILTER_OFF) {
      fpar->filter = false;
    } else if(fcont.cmd == COMMAND::MEAN_FILTER_ON) {
      fpar->mean_filter = true;
    } else if(fcont.cmd == COMMAND::MEAN_FILTER_OFF) {
      fpar->mean_filter = false;
    }
  }
  
  float* in_buf = (float *) inputBuffer;
  float* out_buf = (float *) outputBuffer;

  float current_sample = 0.0;
  short channel = 0;
  for (unsigned int i = 0; i < nBufferFrames * fpar->channels; i++) {
    current_sample = in_buf[i];
    channel = i % fpar->channels;
    //channel = i > nBufferFrames ? 1 : 0;
    
    if(fpar->flippiness > 0.0) {
      current_sample = random_flip(in_buf[i], fpar->flippiness);
    }

    current_sample = fpar->rev.tick(current_sample, channel);
   
    if(fpar->mode == PMODE::MILD){
      // weed out impossible values ... while this kinda takes some of the edge off,
      // it enables us to use some stuff like filters and gain ...
      if (isnan(current_sample)) current_sample = 0;    
      if (current_sample < -1.0) current_sample = -1.0;
      if (current_sample > 1.0) current_sample = 1.0;
    }  
    
    if(fpar->mean_filter){
      fpar->m_fbank->apply(channel, current_sample);
    }
    
    if (fpar->filter) {    
      fpar->fbank->apply(channel, current_sample);      
    } 

    if(fpar->mode == PMODE::MILD){
      current_sample *= fpar->gain;
    } 
  
    out_buf[i] = current_sample;
  }
  
  
  //fpar->rev->processreplace(in_buf, in_buf + nBufferFrames, out_buf, out_buf + nBufferFrames, nBufferFrames);

  return 0;
}

/*
 * INITIALIZE COMMAND LINE OPTIONS !!
 */
// A helper function to simplify the main part (from po example ...)
template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
  return os;
}

// initialize the command line parameter parser
po::options_description init_opts(int ac, char *av[], po::variables_map& vm,
                                  options_container& opts) {

  po::options_description desc("Parameters to use in a creative way");
  desc.add_options()
    ("help", "Display this help!")
    ("input-file", po::value<std::string>(&opts.filename)->default_value(""), "The input file - WAV of FLAC!)")
    ("init-state", po::value<PSTATE>(&opts.initial_state)->default_value(PLAY), "Initial state!")    
    ("init-mode", po::value<PMODE>(&opts.initial_mode)->default_value(MILD), "Initial mode!")
    ("init-gain", po::value<float>(&opts.initial_gain)->default_value(0.5), "Initial gain (default: 0.5)!")
    ("start", po::value<float>(&opts.start)->default_value(0.0), "Starting point within the sample, relative to length!")
    ("end", po::value<float>(&opts.end)->default_value(1.0), "End point within the sample, relative to length!")
    ("sample-repeat", po::value<float>(&opts.sample_repeat)->default_value(1), "Repeat every sample n times!")
    ("buffer-cut", po::value<float>(&opts.buffer_cut)->default_value(2), "Don't fill source output buffer completely!")
    ("offset-mod", po::value<float>(&opts.offset_cut)->default_value(2), "Modify offset increment (chunk size read from buffer)!")
    ("fuzziness", po::value<float>(&opts.fuzziness)->default_value(0.0), "Create fuzziness by removing random samples with a certain probability!")
    ("flippiness", po::value<float>(&opts.flip_prob)->default_value(0.0), "Create different fuzziness by flipping bits a certain probability!")
    ("read-type", po::value<RWTYPES>(&opts.read_type)->default_value(SHORT), "Type used to read from audio file!")
    ("write-type", po::value<RWTYPES>(&opts.write_type)->default_value(SHORT), "Type used to write to audio buffer!")
    ("stream-type", po::value<RWTYPES>(&opts.stream_type)->default_value(SHORT), "Type used for the audio stream!")
    ("channel-offset", po::value<int>(&opts.channel_offset)->default_value(0), "Offset to control channels (esp. useful for mono playback)!")
    ("mean-filter", po::value<int>(&opts.mean_filter_points)->default_value(0), "Apply mean filter to shave the edge off a little!")
    ("mono", po::value<bool>(&opts.mono)->default_value(false), "Mixdown to mono!")
    
    ;
  // ----- end options ... what kind of syntax is this ??

  po::positional_options_description p;
  p.add("input-file", -1);

  po::store(po::command_line_parser(ac, av).options(desc).positional(p).run(),vm);
  po::notify(vm);

  return desc;
}

/*
 * Function to stop audio
 */ 
void stop_audio(RtAudio& source, RtAudio& filter, bool raw) {
  source.stopStream();
  source.closeStream();

  if (!raw){
    filter.stopStream();
    filter.closeStream();
  }
}

template<typename READ_TYPE>
void handle_input(source_params<READ_TYPE>& spar, filter_params& fpar ){
  char input;
  // main loop
  
  while((input = getch()) != 'q'){   
    filter_command_container fcont;
    source_command_container scont;
    switch(input) {    
      // filter control
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
    case '0':
      // stupidly simple, but it seems to work ...
      if (input == '0'){input = 9;}
      else {input -= 49;}
      fpar.fbank->toggle_band(input);
      fpar.fbank->print_bands();
      break;
    case 'f':
      if (fpar.filter) {	
	fcont.cmd = COMMAND::FILTER_OFF;
	std::cout << "filter off" << std::endl;
      } else {
	fcont.cmd = COMMAND::FILTER_ON;
	std::cout << "filter on" << std::endl;
      }
      fpar.cmd_queue->push(fcont);      
      break;
    case 'g':
      if (fpar.mean_filter) {	
	fcont.cmd = COMMAND::MEAN_FILTER_OFF;
	std::cout << "mean (smoothing) filter off" << std::endl;
      } else {
	fcont.cmd = COMMAND::MEAN_FILTER_ON;
	std::cout << "mean (smoothing) filter on" << std::endl;
      }
      fpar.cmd_queue->push(fcont);      
      break;
    // gain control
    case 'd':           
      fcont.cmd = COMMAND::GAIN_CHANGE;
      fcont.gain = fpar.gain + 0.05 >= 1.0 ? 1.0 :  fpar.gain + 0.05;
      std::cout << "gain up, new gain: " << fcont.gain << std::endl;
      fpar.cmd_queue->push(fcont);
      break;
    case 'c':      
      fcont.cmd = COMMAND::GAIN_CHANGE;
      fcont.gain = fpar.gain - 0.05 <= 0.0 ? 0.0 : fpar.gain - 0.05;
      std::cout << "gain down, new gain: " << fcont.gain << std::endl;
      fpar.cmd_queue->push(fcont);
      break;
      // fuzziness control
    case 'a':           
      scont.cmd = COMMAND::FUZZINESS_CHANGE;
      scont.new_fuzziness = spar.fuzziness + 0.05 >= 1.0 ? 1.0 : spar.fuzziness + 0.05;
      std::cout << "fuzziness up, new fuzziness: " << scont.new_fuzziness << std::endl;
      spar.cmd_queue->push(scont);
      break;
    case 'y':      
      scont.cmd = COMMAND::FUZZINESS_CHANGE;
      scont.new_fuzziness = spar.fuzziness - 0.05 <= 0.0 ? 0.0 : spar.fuzziness - 0.05;
      std::cout << "fuzziness down, new fuzziness: " << scont.new_fuzziness << std::endl;
      spar.cmd_queue->push(scont);
      break;
      // samplerate control
    case 's':           
      scont.cmd = COMMAND::SAMPLERATE_CHANGE;
      scont.new_sample_repeat = spar.sample_repeat + 1;
      std::cout << "sample repetition up, now: " << scont.new_sample_repeat << std::endl;
      spar.cmd_queue->push(scont);
      break;
    case 'x':      
      scont.cmd = COMMAND::SAMPLERATE_CHANGE;
      scont.new_sample_repeat = spar.sample_repeat - 1 <= 1 ? 1 : spar.sample_repeat - 1;
      std::cout << "sample repetition down, now: " << scont.new_sample_repeat << std::endl;
      spar.cmd_queue->push(scont);
      break;
      // thread blocking control
    case 'm':      
      scont.cmd = COMMAND::STATE_CHANGE;
      if(spar.state == PLAY){
	scont.new_state = SILENCE;
      } else if (spar.state == LOOP){
	scont.new_state = LOOP_SILENCE;
      } else if (spar.state == SILENCE){
	scont.new_state = PLAY;
      } else if (spar.state == LOOP_SILENCE) {
	scont.new_state = LOOP;
      }
      std::cout << "suspend & mute" << std::endl;
      spar.cmd_queue->push(scont);
      break;    
    case ' ':      
      if(spar.state == PLAY){
	scont.cmd = COMMAND::LOOP_INIT;
	std::cout << "loop from: " << spar.offset << " (~" << (float) spar.offset / spar.fc->samples  << ")" << std::endl;
      } else if (spar.state == LOOP_REC) {
	scont.cmd = COMMAND::LOOP_FINISH;
	std::cout << "loop to: " << spar.offset << " (~" << (float) spar.offset / spar.fc->samples  << ")" << std::endl;
      }
      spar.cmd_queue->push(scont);
      break;
    case 'b':      
      scont.cmd = COMMAND::STATE_CHANGE;
      if(spar.state == PLAY){
	std::cout << "blocking source" << std::endl;
	scont.new_state = BLOCK;
      } else if (spar.state == LOOP){
	std::cout << "blocking source loop" << std::endl;
	scont.new_state = LOOP_BLOCK;
      } else if (spar.state == LOOP_BLOCK) {
	std::cout << "un-blocking source-loop" << std::endl;
	scont.new_state = LOOP;
	cv.notify_all();
      } else if (spar.state == BLOCK){
	std::cout << "un-blocking source" << std::endl;
	scont.new_state = PLAY;
	cv.notify_all();
      }
      spar.cmd_queue->push(scont);
      break;
    default:
      std::cout << "COMMAND NOT ACCEPTED!" << std::endl;
    }
    
  }
}

/*
* The Audio initializing function !
*/
template <typename READ_TYPE, typename WRITE_TYPE>
int handle_audio(options_container& opts) {

  if (opts.filename == "") {
    std::cout << "\nPlease specify input file!\n";
    return EXIT_FAILURE;
  }
  
  source_params<READ_TYPE> spar(opts);

  if (spar.fc->state != file_container<READ_TYPE>::READY) {
    std::cout << "File \"" << opts.filename << "\" is not valid!" << std::endl;
    return EXIT_FAILURE;
  }
  
  spar.fc->print_file_info();
  filter_params fpar(opts, spar.fc->channels);

  std::cout << "Current Parameters:" << std::endl;
  std::cout << "  Read type:     " << rwtypes_strings[opts.read_type] << std::endl;
  std::cout << "  Write type:    " << rwtypes_strings[opts.write_type] << std::endl;
  std::cout << "  Stream type:   " << rwtypes_strings[opts.stream_type] << std::endl;
  std::cout << "  Sample repeat: " << opts.sample_repeat << std::endl;
  std::cout << "  Buffer mod:    " << opts.buffer_cut << std::endl;
  std::cout << "  Offset mod:    " << opts.offset_cut << std::endl << std::endl;
  
  // get current pid
  std::ostringstream str_pid;
  str_pid << ::getpid();

  // initialize audio streams
  RtAudio source;
  RtAudio filter;

  if (filter.getDeviceCount() < 1) {
    std::cout << "\nNo audio devices found!\n";
    return EXIT_FAILURE;
  }
  
  RtAudio::StreamParameters rt_source_parameters;
  
  // in raw mode, all the beatiful glitches will go directly to the DAC ... good luck !
  if(opts.initial_mode != PMODE::RAW){
    RtAudio::StreamParameters rt_filter_out_parameters;
    rt_filter_out_parameters.deviceId = filter.getDefaultOutputDevice();
    rt_filter_out_parameters.nChannels = spar.fc->channels;
    rt_filter_out_parameters.firstChannel = opts.channel_offset;

    RtAudio::StreamParameters rt_filter_in_parameters;
    rt_filter_in_parameters.deviceId = filter.getDefaultOutputDevice();
    rt_filter_in_parameters.nChannels = spar.fc->channels;
    rt_filter_in_parameters.firstChannel = opts.channel_offset;

    RtAudio::StreamOptions rt_filter_opts;
    // rt_filter_opts.flags = RTAUDIO_NONINTERLEAVED;
    rt_filter_opts.streamName = "akita-filter-";
    rt_filter_opts.streamName.append(str_pid.str());
    rt_filter_opts.autoConnectInput = false;

    // initialize filter client
    try {
      filter.openStream(&rt_filter_out_parameters, &rt_filter_in_parameters, RTAUDIO_FLOAT32, opts.samplerate,
			&opts.buffer_frames, &filter_callback,
			(void *)&fpar, &rt_filter_opts);
      filter.startStream();
    } catch (RtAudioError &e) {
      e.printMessage();
      return EXIT_FAILURE;
    }
    
    //find filter thread
    for(unsigned int i = 0; i < source.getDeviceCount(); i++){
      if (source.getDeviceInfo(i).name == rt_filter_opts.streamName){
	rt_source_parameters.deviceId = i;
      }
    }
  } else {
    rt_source_parameters.deviceId = source.getDefaultOutputDevice();;
  }

  rt_source_parameters.nChannels = spar.fc->channels;
  rt_source_parameters.firstChannel = 0;

  RtAudio::StreamOptions source_opts;
  source_opts.streamName = "akita-noise-source-";
  source_opts.streamName.append(str_pid.str());

  // initialize noise source
  try {
    source.openStream(&rt_source_parameters, NULL, rwtypes_rtaudio[opts.stream_type], opts.samplerate,
		       &opts.buffer_frames, &source_callback<READ_TYPE, WRITE_TYPE>,
		       (void *)&spar, &source_opts);
    
    source.startStream();
  } catch (RtAudioError &e) {
    e.printMessage();
    return EXIT_FAILURE;
  }

  std::cout << "Playing ... press q to quit" << std::endl;
  
  handle_input(spar, fpar);
  
  try {
    stop_audio(source, filter, (opts.initial_mode == PMODE::RAW));
  } catch(RtAudioError &e){
    e.printMessage();
    return EXIT_FAILURE;
  }

  std::cout << "\nQuitting, bye!.\n";
  
  return EXIT_SUCCESS;
}

typedef std::pair<RWTYPES, RWTYPES> init_key;

typedef std::function<int (options_container&)> audio_handler;

std::map<init_key, audio_handler> handlers {
  { init_key(SHORT, SHORT), handle_audio<int16_t, int16_t> },
  { init_key(SHORT, FLOAT), handle_audio<int16_t, float> },
  { init_key(SHORT, DOUBLE), handle_audio<int16_t, double> },
  { init_key(SHORT, UCHAR), handle_audio<int16_t, uint8_t> },
  { init_key(FLOAT, SHORT), handle_audio<float, int16_t> },
  { init_key(FLOAT, FLOAT), handle_audio<float, float> },
  { init_key(FLOAT, DOUBLE), handle_audio<float, double> },
  { init_key(FLOAT, UCHAR), handle_audio<float, uint8_t> },
  { init_key(DOUBLE, SHORT), handle_audio<double, int16_t> },
  { init_key(DOUBLE, FLOAT), handle_audio<double, float> },
  { init_key(DOUBLE, DOUBLE), handle_audio<double, double> },
  { init_key(DOUBLE, UCHAR), handle_audio<double, uint8_t> },
};

/*
 * MAIN ROUTINE!!
 */
int main(int ac, char *av[]) {
  // Salutations!
  std::cout << "\n~~ akita - create noise abusing low-level audio parameters! ~~\n" << std::endl;

  options_container opts;
  
   // initialize command line options
  po::variables_map vm;
  po::options_description desc = init_opts(ac, av, vm, opts);

  // help output
  if (vm.count("help")) {
    std::cout << "usage: akita <file> [options] \n" << std::endl;
    std::cout << desc;
    std::cout << "\nRead/Write/Stream types options:" << std::endl;
    std::cout << "  uchar      8 Bit, unsigned char (not possible as read type!)" << std::endl;
    std::cout << "  short      16 Bit, short (default)" << std::endl;
    std::cout << "  float      32 Bit, float" << std::endl;
    std::cout << "  double     64 Bit, double" << std::endl;
    return EXIT_SUCCESS;
  }
   
  //from here on, things must be parametrized
  return handlers[init_key(opts.read_type, opts.write_type)](opts);;
}
