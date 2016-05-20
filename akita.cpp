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
  enum COMMAND { MODE_CHANGE, STATE_CHANGE, GAIN_CHANGE, LOOP_INIT, LOOP_FINISH, LOOP_RELEASE };
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
  
  // format glitch parameters
  RWTYPES read_type;
  RWTYPES write_type;
  RWTYPES stream_type;

  // buffer glitch parameters
  float buffer_cut;
  float offset_cut;
  float sample_repeat;

  float initial_gain;
};

// file container to load and store samples in buffer
template <typename READ_TYPE>
struct file_container {

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
    
    // libsndfile soundfilehandle ... 
    SndfileHandle file = SndfileHandle(fname.c_str());

    // transfer some info
    samplerate = file.samplerate();
    channels = file.channels();
    frames = file.frames();
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
  }

  ~file_container() { delete[] file_buffer; }

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
};

struct filter_command_container {
  COMMAND::COMMAND cmd;
  
  // parameters that could be subject to change
  PMODE new_mode;
  float gain;
  bool filter;
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

  source_params(options_container& opts) {
    state = opts.initial_state;
    buffer_cut = opts.buffer_cut;
    offset_cut = opts.offset_cut;
    sample_repeat = opts.sample_repeat;
    
    fc = new file_container<READ_TYPE>(opts.filename);
    cmd_queue = new lfree::spsc_queue<source_command_container>(10);
  }

  ~source_params(){
    delete fc;
    delete cmd_queue;
  }
};

// params for the filter 
struct filter_params {
  lfree::spsc_queue<filter_command_container>* cmd_queue;
  PMODE mode;
  float gain;
  bool filter = false;

  filter_params (options_container& opts) {
    mode = opts.initial_mode;
    gain = opts.initial_gain;

    cmd_queue = new lfree::spsc_queue<filter_command_container>(10);
  }

  ~filter_params () {
    delete cmd_queue;
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
  source_command_container cont;
  while(cmd_queue->pop(&cont)){
    if(cont.cmd == COMMAND::STATE_CHANGE){
      spar->state = cont.new_state;
    } else if(cont.cmd == COMMAND::LOOP_INIT) {
      spar->loop_start = spar->offset;
      spar->state = LOOP_REC;
    } else if(cont.cmd == COMMAND::LOOP_FINISH) {
      spar->loop_end = spar->offset;
      spar->state = LOOP;
    }
  }

  // block source thread ... might have an interesting effect ...
  if(spar->state == BLOCK){
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
  for (int i = 0; i < nBufferFrames * spar->buffer_cut; i++) {
    for (float j = 0; j < spar->sample_repeat; j += 1.0) {
      if (spar->offset + i > fc->samples) {
        spar->offset = 0;
      }
      *out_buf++ = fc->file_buffer[spar->offset + i];          
    }
  }

  // increament read offest
  spar->offset += nBufferFrames * spar->offset_cut;
  if((spar->state == LOOP) && (spar->offset >= spar->loop_end)){
    spar->offset = spar->loop_start;
  }
 
  return 0;
}


// the callback function for the filter thread
int filter_callback(void *outputBuffer, void *inputBuffer,
                 unsigned int nBufferFrames, double streamTime,
	   RtAudioStreamStatus status, void *userData) {

  float* in_buf = (float *) inputBuffer;
  float* out_buf = (float *) outputBuffer;
  
  for (unsigned int i = 0; i < nBufferFrames * 2; i++) {

    float current_sample = in_buf[i];

    if (isnan(current_sample)) current_sample = 0;    
    if (current_sample < -1.0) current_sample = -1.0;
    if (current_sample > 1.0) current_sample = 1.0;
        
    out_buf[i] = current_sample;
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
    ("init-mode", po::value<PMODE>(&opts.initial_mode)->default_value(MILD), "Initial state!")
    ("init-gain", po::value<float>(&opts.initial_gain)->default_value(0.5), "Initial gain (default: 0.5)!")
    ("sample-repeat", po::value<float>(&opts.sample_repeat)->default_value(1), "Repeat every sample n times!")
    ("buffer-mod", po::value<float>(&opts.buffer_cut)->default_value(2), "Don't fill source output buffer completely!")
    ("offset-mod", po::value<float>(&opts.offset_cut)->default_value(2), "Modify offset increment (chunk size read from buffer)!")
    ("read-type", po::value<RWTYPES>(&opts.read_type)->default_value(SHORT), "Type used to read from audio file!")
    ("write-type", po::value<RWTYPES>(&opts.write_type)->default_value(SHORT), "Type used to write to audio buffer!")
    ("stream-type", po::value<RWTYPES>(&opts.stream_type)->default_value(SHORT), "Type used for the audio stream!")
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
void stop_audio(RtAudio& source, RtAudio& filter) {
  
  source.stopStream();
  source.closeStream();

  filter.stopStream();
  filter.closeStream();
}

template<typename READ_TYPE>
void handle_input(source_params<READ_TYPE> spar, filter_params fpar ){
  char input;
  // main loop
  
  while((input = getch()) != 'q'){   
    source_command_container cont;
    switch(input) {
    case 'd':
      std::cout << "gain up" << std::endl;
      //cont.cmd = COMMAND::GAIN_CHANGE;      
      break;
    case 'c':
      std::cout << "gain down" << std::endl;
      //   cont.cmd = COMMAND::GAIN_CHANGE;      
      break;
    case 's':
      cont.cmd = COMMAND::STATE_CHANGE;
      if(spar.state == PLAY){
	cont.new_state = SILENCE;
      } else if (spar.state == LOOP){
	cont.new_state = LOOP_SILENCE;
      } else if (spar.state == SILENCE){
	cont.new_state = PLAY;
      } else if (spar.state == LOOP_SILENCE) {
	cont.new_state = LOOP;
      }
      break;    
    case ' ':
      if(spar.state == PLAY){
	cont.cmd = COMMAND::LOOP_INIT;
	std::cout << "loop from: " << spar.offset << std::endl;
      } else if (spar.state == LOOP_REC) {
	cont.cmd = COMMAND::LOOP_FINISH;
	std::cout << "loop to: " << spar.offset << std::endl;
      }
      break;
    case 'b':
      cont.cmd = COMMAND::STATE_CHANGE;
      if(spar.state == PLAY){
	cont.new_state = BLOCK;
      } else if (spar.state == LOOP){
	cont.new_state = LOOP_BLOCK;
      } else if (spar.state == LOOP_BLOCK) {
	cont.new_state = LOOP;
	cv.notify_all();
      } else if (spar.state == BLOCK){
	cont.new_state = PLAY;
	cv.notify_all();
      }
      break;
    default:
      std::cout << "COMMAND NOT ACCEPTED!" << std::endl;
    }
    spar.cmd_queue->push(cont);
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
  
  RtAudio source;
  RtAudio filter;

  if (filter.getDeviceCount() < 1) {
    std::cout << "\nNo audio devices found!\n";
    return EXIT_FAILURE;
  }
  
  source_params<READ_TYPE> spar(opts);
  spar.fc->print_file_info();
  filter_params fpar(opts);
  
  // get current pid
  std::ostringstream str_pid;
  str_pid << ::getpid();
  
  RtAudio::StreamParameters rt_filter_out_parameters;
  rt_filter_out_parameters.deviceId = filter.getDefaultOutputDevice();
  rt_filter_out_parameters.nChannels = 2;
  rt_filter_out_parameters.firstChannel = 0;

  RtAudio::StreamParameters rt_filter_in_parameters;
  rt_filter_in_parameters.deviceId = filter.getDefaultOutputDevice();
  rt_filter_in_parameters.nChannels = 2;
  rt_filter_in_parameters.firstChannel = 0;

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

  RtAudio::StreamParameters rt_source_parameters;

  //find filter thread
  for(unsigned int i = 0; i < source.getDeviceCount(); i++){
    if (source.getDeviceInfo(i).name == rt_filter_opts.streamName){
      rt_source_parameters.deviceId = i;
    }
  }
  
  rt_source_parameters.nChannels = 2;
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

  std::cout << "\nPlaying ... press q to quit.\n";
  
  handle_input(spar, fpar);
  
  try {
    stop_audio(source, filter);
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
    
  std::cout << "Current Parameters:" << std::endl;
  std::cout << "  Read type:     " << rwtypes_strings[opts.read_type] << std::endl;
  std::cout << "  Write type:    " << rwtypes_strings[opts.write_type] << std::endl;
  std::cout << "  Stream type:   " << rwtypes_strings[opts.stream_type] << std::endl;
  std::cout << "  Sample repeat: " << opts.sample_repeat << std::endl;
  std::cout << "  Buffer mod:    " << opts.buffer_cut << std::endl;
  std::cout << "  Offset mod:    " << opts.offset_cut << std::endl << std::endl;

  //from here on, things must be parametrized
  int exit_code = handlers[init_key(opts.read_type, opts.write_type)](opts);

  std::cout << exit_code << std::endl;
  
  return exit_code;
}
