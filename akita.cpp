#include <iostream>
#include <algorithm>
#include <iterator>
#include <sndfile.hh>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include "getch.h"
#include "RtAudio.h"
#include <cmath>
#include <mutex>              
#include <condition_variable>
#include <string>
#include <functional>
#include <climits>

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
  enum COMMAND { MODE_CHANGE, STATE_CHANGE, SAMPLERATE_CHANGE, FUZZINESS_CHANGE, GAIN_CHANGE, FILTER_ON, FILTER_OFF, LOOP_INIT, LOOP_FINISH, LOOP_RELEASE };
}

/*
 * Filter and filterbank
 */
struct canonical_sos_filter {
  enum FMODE {HP, BP, LP, NOTCH, AP};

  float a1, a2;
  float b0, b1, b2;
  float k;
  float q;

  float del1, del2; 
  
  canonical_sos_filter(){
    update(5000, 2, 44100, LP);
  }

  canonical_sos_filter(double frequency, double q, int samplerate, FMODE mode){
    update(frequency, q, samplerate, mode);
  }

  void update(double frequency, double q, int samplerate, FMODE mode) {    
    del1 = 0;
    del2 = 0;
    k = tanh( (M_PI*frequency) / samplerate);
    a1 = (2.0 * q * (pow(k,2) - 1)) / ((pow(k,2) * q) + k + q);
    a2 = ((pow(k,2) * q) - k + q) / ((pow(k,2) * q) + k + q);
    if (mode == FMODE::LP){
      b0 = (pow(k,2) * q) / ((pow(k,2) * q) + k + q);
      b1 = (2.0 * pow(k,2) * q) / ((pow(k,2) * q) + k + q);
      b2 = b0;      
    } else if (mode == FMODE::HP){
      b0 = q / ((pow(k,2) * q) + k + q);
      b1 = -1.0 * ((2.0 * q) / ((pow(k,2) * q) + k + q));
      b2 = b0;      
    } else if (mode == FMODE::BP){
      b0 = k / ((pow(k,2) * q) + k + q);
      b1 = 0;
      b2 = -1.0 * b0;      
    } else if (mode == FMODE::NOTCH){
      b0 = (q * (1.0 + pow(k,2))) / ((pow(k,2) * q) + k + q);
      b1 = (2 * q * (pow(k,2) - 1)) / ((pow(k,2) * q) + k + q);
      b2 = b0;      
    } else if (mode == FMODE::AP){
      b0 = ((pow(k,2) * q) - k + q) / ((pow(k,2) * q) + k + q);
      b1 = (2 * q * (pow(k,2) - 1)) / ((pow(k,2) * q) + k + q);
      b2 = 1.0;      
    }
  }

  double calculate(float sample){
    float intermediate = sample + ((-1.0 * a1) * del1) + ((-1.0 * a2) * del2);
    float out = (b0 * intermediate) + (b1 * del1) + (b2 * del2);
    del2 = del1;
    del1 = intermediate;
    return out;
  }

  void process(float& sample){
    sample = calculate(sample);          
  }
};

// a simple filterbank consisting of several state-variable filters
struct filterbank {
  filterbank(int channels, int samplerate, double lowcut, double hicut, int bands) {
    this->bands = bands;
    this->channels = channels;
    
    fbank = new canonical_sos_filter[channels * bands];

    for(int ch = 0; ch < channels; ch++) {
      fbank[ch * bands].update(lowcut, 1.5, samplerate, canonical_sos_filter::HP);
      for(int b = 1; b < bands-1; b++) {
	fbank[(ch * bands) + b].update(lowcut + (b * ((hicut - lowcut) / bands)) , 1.5, samplerate, canonical_sos_filter::NOTCH);
      }
      fbank[(ch * bands) + bands-1].update(hicut, 1.5, samplerate, canonical_sos_filter::LP);
    }

    fmask = new bool[bands];
    for(int b = 0; b < bands; b++) {
      fmask[b] = false;
    }
  }

  ~filterbank() {
    delete [] fbank;
    delete [] fmask;
  }

  int bands;
  int channels;
  canonical_sos_filter* fbank;
  bool* fmask;

  float apply(int channel, float& sample) {    
    for(int b = 0; b < bands; b++){
      if(fmask[b]){	
	fbank[(channel * bands) + b].process(sample);
      }
    }
    return sample;
  }

  void toggle_band(int band) {
    if(fmask[band]) {
      fmask[band] = false;
    } else {
      fmask[band] = true;
    } 
   }

  void print_bands(){
    std::cout << "Filter bands: ";
    for(int i = 0; i < bands; i++){
      std::cout << "[" << fmask[i] << "] ";
    }
    std::cout << std::endl;
  }
};

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
  
  // format glitch parameters
  RWTYPES read_type;
  RWTYPES write_type;
  RWTYPES stream_type;

  // buffer glitch parameters
  float buffer_cut;
  float offset_cut;
  float sample_repeat;

  float initial_gain;

  int channel_offset;
  
  // kill samples to make it all fuzzy !
  float fuzziness;
};

struct file_container_exception{
  std::string msg;
  file_container_exception(std::string message){
    msg = message;
  }
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
  short channels;
  int samplerate;

  // methods
  file_container(std::string fname) {
    name = fname;      
  }

  // load file to memory
  void load(){
    // libsndfile soundfilehandle ... 
    SndfileHandle file = SndfileHandle(name.c_str());

    channels = file.channels();
    // transfer some info
    samplerate = file.samplerate();

    frames = file.frames();

    // can't work on an empty file !
    if(frames <= 0){
      throw file_container_exception("Soundfile \"" + name + "\" invalid !");
    }
    
    samples = frames * channels;
  
    file_buffer = new READ_TYPE[samples];

    int frame_chunksize = CHUNKSIZE * channels;

    READ_TYPE *chunk_buffer = new READ_TYPE[frame_chunksize];

    // Read file into buffer -- could this be made faster ?
    long int read_offset = 0;
    while (file.readf(chunk_buffer, CHUNKSIZE) == CHUNKSIZE) {
      for (int i = 0; i < frame_chunksize; i++) {
	file_buffer[read_offset + i] = chunk_buffer[i];
      }
      read_offset += frame_chunksize;
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
    
    fc = new file_container<READ_TYPE>(opts.filename);
    buffer_cut = fc->channels;
    offset_cut = fc->channels;
    opts.buffer_cut = buffer_cut;
    opts.offset_cut = offset_cut;
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
  PMODE mode;
  int channels;
  float gain;
  bool filter = false;
  
  filter_params (options_container& opts, int channels) {
    mode = opts.initial_mode;
    gain = opts.initial_gain;
    this->channels = channels;
    
    cmd_queue = new lfree::spsc_queue<filter_command_container>(10);
    fbank = new filterbank(channels, opts.samplerate, 100, 8000, 10);
  }

  ~filter_params () {
    //std::cout << "cleaning filter params" << std::endl;
    delete cmd_queue;
    delete fbank;
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
    if (spar->offset + i > fc->samples) {
      spar->offset = 0;
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
    }
  }
  
  float* in_buf = (float *) inputBuffer;
  float* out_buf = (float *) outputBuffer;

  if(fpar->mode == PMODE::MILD){
    for (unsigned int i = 0; i < nBufferFrames * fpar->channels; i++) {

      float current_sample = in_buf[i];

      // weed out impossible values ... while this kinda takes some of the edge off,
      // it enables us to use some stuff like filters and gain ...
      if (isnan(current_sample)) current_sample = 0;    
      if (current_sample < -1.0) current_sample = -1.0;
      if (current_sample > 1.0) current_sample = 1.0;
      if (fpar->filter) {
	current_sample *= fpar->gain;
	out_buf[i] = fpar->fbank->apply(i % fpar->channels, current_sample);	  	
      } else {
	out_buf[i] = current_sample * fpar->gain;	  	
      }
    }
  } else {
    // in pure prox mode, all kinds of manipulations hardly make sense
    for (unsigned int i = 0; i < nBufferFrames * fpar->channels; i++) {
      out_buf[i] = in_buf[i];
    }
  }

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
    ("sample-repeat", po::value<float>(&opts.sample_repeat)->default_value(1), "Repeat every sample n times!")
    ("buffer-cut", po::value<float>(&opts.buffer_cut)->default_value(2), "Don't fill source output buffer completely!")
    ("offset-mod", po::value<float>(&opts.offset_cut)->default_value(2), "Modify offset increment (chunk size read from buffer)!")
    ("fuzziness", po::value<float>(&opts.fuzziness)->default_value(0.0), "Create fuzziness by removing random samples with a certain probability!")
    ("read-type", po::value<RWTYPES>(&opts.read_type)->default_value(SHORT), "Type used to read from audio file!")
    ("write-type", po::value<RWTYPES>(&opts.write_type)->default_value(SHORT), "Type used to write to audio buffer!")
    ("stream-type", po::value<RWTYPES>(&opts.stream_type)->default_value(SHORT), "Type used for the audio stream!")
    ("channel-offset", po::value<int>(&opts.channel_offset)->default_value(0), "Offset to control channels (esp. useful for mono playback)!")
    
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

  //load audio file 
  try{
    spar.fc->load();
  } catch (file_container_exception& e) {
    std::cout << e.msg << std::endl;
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
