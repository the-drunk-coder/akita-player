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

namespace po = boost::program_options;
namespace lfree = boost::lockfree;

#define CHUNKSIZE 1024

// format enum
enum RWTYPES { UCHAR, SHORT, FLOAT, DOUBLE };

namespace PSTATE{
  enum PSTATE { PLAY, LOOP_REC, LOOP, LOOP_SILENCE, LOOP_BLOCK, SILENCE, BLOCK };
}

namespace COMMAND {
  enum COMMAND { STATE_CHANGE, GAIN_CHANGE, LOOP_INIT, LOOP_FINISH };
}

// map enable nicer type output ...
std::map <RWTYPES, std::string> rwtypes_strings {
  { UCHAR, "uchar (8 Bit)" },
  { SHORT, "short (16 Bit)" },
  { FLOAT, "float (32 Bit)" },
  { DOUBLE, "double (64 Bit)" },
};

// parameters struct
struct params_to_abuse {
  float buffer_cut;
  float offset_cut;
  float sample_repeat;
  RWTYPES read_type;
  RWTYPES write_type;
  RWTYPES stream_type;
};

struct command_container {
  COMMAND::COMMAND cmd;
  PSTATE::PSTATE new_state;
  float new_gain;
};

// naive state variable filter, modeled after DAFx Chapter 2
struct state_variable_filter {
  state_variable_filter(double frequency, double q, int samplerate){
    q1 = 1.0 / q;
    f1 = 2 * sin(M_PI * frequency / samplerate);
    del_lp = 0.0;
    del_bp = 0.0;
    hp = 0.0;
    lp = 0.0;
    bp = 0.0;
  }

  //double frequency;
  //double q;
  double q1;
  double f1;

  //delays
  double del_lp;
  double del_bp;

  // current
  double hp;
  double bp;
  double lp;

  void calculate(double sample){
    hp = sample - del_lp - q1 * del_bp;
    bp = f1 * hp + del_bp;
    lp = f1 * bp + del_lp;
    del_bp = bp;
    del_lp = lp;
  }
};


struct play_params {
  // loop params
  lfree::spsc_queue<command_container>* cmd_queue;
  PSTATE::PSTATE state;
  state_variable_filter* filter_l;
  state_variable_filter* filter_r; 
  double gain;
  long int loop_start;
  long int loop_end;
  long int offset;
};

// file container to communicate with callback function
template <typename READ_TYPE>
class file_container {
public:
  file_container(SndfileHandle file) {
    samples = file.frames() * file.channels();
    file_buffer = new READ_TYPE[samples];
    channels = file.channels();
  }

  ~file_container() { delete[] file_buffer; }

  READ_TYPE *file_buffer;
  long int samples;
  short channels;

  // the params to make things weird!
  params_to_abuse *params; 

  play_params *pp;
};

// this is bad ... 
std::mutex mtx;
std::condition_variable cv;

// the parameterized callback function ...
template <typename READ_TYPE, typename WRITE_TYPE>
int abusive_play(void *outputBuffer, void *inputBuffer,
                 unsigned int nBufferFrames, double streamTime,
                 RtAudioStreamStatus status, void *userData) {

  // get the parameter container from the user data ...
  file_container<READ_TYPE> *fc = reinterpret_cast<file_container<READ_TYPE> *>(userData);
  params_to_abuse *params = fc->params;

  play_params *pp = fc->pp;
  lfree::spsc_queue<command_container> *cmd_queue = pp->cmd_queue;
  WRITE_TYPE *out_buf = (WRITE_TYPE *)outputBuffer;

  // handle commands
  command_container cont;
  while(cmd_queue->pop(&cont)){
    //std::cout << "rec cmd" << std::endl;
    if(cont.cmd == COMMAND::STATE_CHANGE){
      pp->state = cont.new_state;
    } else if(cont.cmd == COMMAND::GAIN_CHANGE){
      //std::cout << "gain" << std::endl;
      pp->gain = cont.new_gain;
      if(pp->gain > 1.0){
	pp->gain = 1.0;
      } else if (pp->gain < 0.0)  {
	pp->gain = 0.0;
      }
    } else if(cont.cmd == COMMAND::LOOP_INIT) {
      pp->loop_start = pp->offset;
      pp->state = PSTATE::LOOP_REC;
    } else if(cont.cmd == COMMAND::LOOP_FINISH) {
      pp->loop_end = pp->offset;
      pp->state = PSTATE::LOOP;
    }
  }

  if(pp->state == PSTATE::BLOCK){
    std::unique_lock<std::mutex> lck(mtx);
    cv.wait(lck);
    return 0;
  }

  if (status) { std::cout << "Stream underflow detected!" << std::endl; }

  // SILENCE ! I KILL YOU !!
  if (pp->state == PSTATE::SILENCE || pp->state == PSTATE::LOOP_SILENCE) {
    for (int i = 0; i < nBufferFrames * params->buffer_cut; i++) {
      *out_buf++ = 0;
    }
    // in this case, offset won't be modified
    return 0;
  }

  short filter_flag = 0;
  READ_TYPE out_sample;
  // transfer samples from file buffer to output !
  for (int i = 0; i < nBufferFrames * params->buffer_cut; i++) {
    for (float j = 0; j < params->sample_repeat; j += 1.0) {
      if (pp->offset + i > fc->samples) {
        pp->offset = 0;
      }
      *out_buf++ = fc->file_buffer[pp->offset + i];
      /*
      if(filter_flag == 0) {
	pp->filter_l->calculate((double) out_sample);
	*out_buf++ = (WRITE_TYPE) (pp->filter_l->lp);
	filter_flag = 1;
      } else {
	pp->filter_r->calculate((double) out_sample);
	*out_buf++ = (WRITE_TYPE) (pp->filter_r->lp);
	filter_flag = 0;
	} 
      */    
    }
  }

  // increament read offest
  pp->offset += nBufferFrames * params->offset_cut;
  if((pp->state == PSTATE::LOOP) && (pp->offset >= pp->loop_end)){
    pp->offset = pp->loop_start;
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

// custom stream to extract enum from command line
std::istream &operator>>(std::istream &in, RWTYPES &rwtype) {
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

// initialize the command line parameter parser
po::options_description init_opts(int ac, char *av[], po::variables_map *vm,
                                  params_to_abuse *params) {

  po::options_description desc("Parameters to use in a creative way");
  desc.add_options()
    ("help", "Display this help!")
    ("input-file", po::value<std::string>(), "The input file - WAV of FLAC!")
    ("sample-repeat", po::value<float>(), "Repeat every sample n times!")
    ("buffer-mod", po::value<float>(), "Don't fill output buffer completely!")
    ("offset-mod", po::value<float>(), "Modify offset increment (chunk size read from buffer)!")
    ("read-type", po::value<RWTYPES>(&(*params).read_type)->default_value(SHORT), "Type used to read from audio file!")
    ("write-type", po::value<RWTYPES>(&(*params).write_type)->default_value(SHORT), "Type used to write to audio buffer!")
    ("stream-type", po::value<RWTYPES>(&(*params).stream_type)->default_value(SHORT), "Type used for the audio stream!")
    ;
  // ----- end options ... what kind of syntax is this ??

  po::positional_options_description p;
  p.add("input-file", -1);

  po::store(po::command_line_parser(ac, av).options(desc).positional(p).run(),*vm);
  po::notify(*vm);

  return desc;
}

template <typename READ_TYPE, typename WRITE_TYPE>
int load_file_and_play(SndfileHandle file, params_to_abuse *params, play_params *pp, RtAudio *dac) {

  file_container<READ_TYPE> *fc = new file_container<READ_TYPE>(file);
  fc->params = params;
  fc->pp = pp;

  int frame_chunksize = CHUNKSIZE * file.channels();
  READ_TYPE *chunk_buffer = new READ_TYPE[frame_chunksize];

  // Read file into buffer
  long int read_offset = 0;
  while (file.readf(chunk_buffer, CHUNKSIZE) == CHUNKSIZE) {
    for (int i = 0; i < frame_chunksize; i++) {
      fc->file_buffer[read_offset + i] = chunk_buffer[i];
    }
    read_offset += frame_chunksize;
  }
  // not needed any longer !
  delete[] chunk_buffer;

  if (dac->getDeviceCount() < 1) {
    std::cout << "\nNo audio devices found!\n";
    return 0;
  }

  RtAudio::StreamParameters parameters;
  parameters.deviceId = dac->getDefaultOutputDevice();
  parameters.nChannels = 2;
  parameters.firstChannel = 0;
  unsigned int sampleRate = 44100;
  unsigned int bufferFrames = 1024; // 256 sample frames

  try {
    if (params->stream_type == UCHAR) {
      dac->openStream(&parameters, NULL, RTAUDIO_SINT8, sampleRate,
		      &bufferFrames, &abusive_play<READ_TYPE, WRITE_TYPE>,
		      (void *)fc);
    } else if (params->stream_type == FLOAT) {
      dac->openStream(&parameters, NULL, RTAUDIO_FLOAT32, sampleRate,
		      &bufferFrames, &abusive_play<READ_TYPE, WRITE_TYPE>,
		      (void *)fc);
    } else if (params->stream_type == DOUBLE) {
      dac->openStream(&parameters, NULL, RTAUDIO_FLOAT64, sampleRate,
		      &bufferFrames, &abusive_play<READ_TYPE, WRITE_TYPE>,
		      (void *)fc);
    } else {
      // this is probably the most commom option
      dac->openStream(&parameters, NULL, RTAUDIO_SINT16, sampleRate,
		      &bufferFrames, &abusive_play<READ_TYPE, WRITE_TYPE>,
		      (void *)fc);
    }
    dac->startStream();
  } catch (RtAudioError &e) {
    e.printMessage();
    return 0;
  }

  return 0;
}

/*
 * MAIN ROUTINE!!
 */
int main(int ac, char *av[]) {
  // Salutations!
  std::cout << "\n~~ akita - create noise abusing low-level audio parameters! ~~\n" << std::endl;

  // get params to abuse ...
  params_to_abuse params;

  // initialize command line options
  po::variables_map vm;
  po::options_description desc = init_opts(ac, av, &vm, &params);

  // help output
  if (vm.count("help")) {
    std::cout << "usage: akita <file> [options] \n" << std::endl;
    std::cout << desc;
    std::cout << "\nRead/Write/Stream types options:" << std::endl;
    std::cout << "  uchar      8 Bit, unsigned char (not possible as read type!)" << std::endl;
    std::cout << "  short      16 Bit, short (default)" << std::endl;
    std::cout << "  float      32 Bit, float" << std::endl;
    std::cout << "  double     64 Bit, double" << std::endl;
    return 0;
  }

  // get filename to use ...
  const char *fname = "";
  if (vm.count("input-file")) {
    fname = vm["input-file"].as<std::string>().c_str();
  } else {
    std::cout << "Please specify input file !" << std::endl;
    return 0;
  }

  SndfileHandle file = SndfileHandle(fname);

  params.sample_repeat = 1.0;
  if (vm.count("sample-repeat")) {
    params.sample_repeat = vm["sample-repeat"].as<float>();
  }

  params.buffer_cut = (float)file.channels();
  if (vm.count("buffer-mod")) {
    params.buffer_cut = vm["buffer-mod"].as<float>();
  }

  params.offset_cut = (float)file.channels();
  if (vm.count("offset-mod")) {
    params.offset_cut = vm["offset-mod"].as<float>();
  }

  // display some file info ...
  std::cout << "Input file: " << fname << ", " << file.channels() << "ch, "
            << file.samplerate() << "Hz, " << file.frames() << " frames.\n"
            << std::endl;

  std::cout << "Current Parameters:" << std::endl;
  std::cout << "  Read type:     " << rwtypes_strings[params.read_type] << std::endl;
  std::cout << "  Write type:    " << rwtypes_strings[params.write_type] << std::endl;
  std::cout << "  Stream type:   " << rwtypes_strings[params.stream_type] << std::endl;
  std::cout << "  Sample repeat: " << params.sample_repeat << std::endl;
  std::cout << "  Buffer mod:    " << params.buffer_cut << std::endl;
  std::cout << "  Offset mod:    " << params.offset_cut << std::endl;

  lfree::spsc_queue<command_container> cmd_queue(10);

  // play params, for looping etc ...
  play_params pp;
  state_variable_filter filter_l(3000, 7, file.samplerate());
  state_variable_filter filter_r(3000, 7, file.samplerate());
  
  pp.gain = 0.001;
  pp.filter_l = &filter_l;
  pp.filter_r = &filter_r;
  pp.state = PSTATE::PLAY;
  pp.cmd_queue = &cmd_queue;
  pp.loop_start = 0;
  pp.loop_end = 0;
  pp.offset = 0;

  // ok, here it get's a little awkward ... a dynamically-typed language would
  // come
  // in handy here ...
  int exit_code = 0;

  RtAudio dac;
  if (params.read_type == FLOAT) {
    if (params.write_type == UCHAR) {
      exit_code = load_file_and_play<float_t, int8_t>(file, &params, &pp, &dac);
    } else if (params.write_type == FLOAT) {
      exit_code = load_file_and_play<float_t, float_t>(file, &params, &pp, &dac);
    } else if (params.write_type == DOUBLE) {
      exit_code = load_file_and_play<float_t, double_t>(file, &params, &pp, &dac);
    } else {
      exit_code = load_file_and_play<float_t, int16_t>(file, &params, &pp, &dac);
    }
  } else if (params.read_type == DOUBLE) {
    if (params.write_type == UCHAR) {
      exit_code = load_file_and_play<double_t, int8_t>(file, &params, &pp, &dac);
    } else if (params.write_type == FLOAT) {
      exit_code = load_file_and_play<double_t, float_t>(file, &params, &pp, &dac);
    } else if (params.write_type == DOUBLE) {
      exit_code = load_file_and_play<double_t, double_t>(file, &params, &pp, &dac);
    } else {
      exit_code = load_file_and_play<double_t, int16_t>(file, &params, &pp, &dac);
    }
  } else {
    if (params.write_type == UCHAR) {
      exit_code = load_file_and_play<int16_t, int8_t>(file, &params, &pp, &dac);
    } else if (params.write_type == FLOAT) {
      exit_code = load_file_and_play<int16_t, float_t>(file, &params, &pp, &dac);
    } else if (params.write_type == DOUBLE) {
      exit_code = load_file_and_play<int16_t, double_t>(file, &params, &pp, &dac);
    } else {
      exit_code = load_file_and_play<int16_t, int16_t>(file, &params, &pp, &dac);
    }
  }

  char input;
  std::cout << "\nPlaying ... press q to quit.\n";
  // main loop
  while((input = getch()) != 'q'){   
    command_container cont;
    switch(input) {
    case 'd':
      std::cout << "gain up" << std::endl;
      cont.cmd = COMMAND::GAIN_CHANGE;
      cont.new_gain = pp.gain + 0.001;      
      break;
    case 'c':
      std::cout << "gain down" << std::endl;
      cont.cmd = COMMAND::GAIN_CHANGE;      
      cont.new_gain = pp.gain - 0.001;
      break;
    case 's':
      cont.cmd = COMMAND::STATE_CHANGE;
      if(pp.state == PSTATE::PLAY){
	cont.new_state = PSTATE::SILENCE;
      } else if (pp.state == PSTATE::LOOP){
	cont.new_state = PSTATE::LOOP_SILENCE;
      } else if (pp.state == PSTATE::SILENCE){
	cont.new_state = PSTATE::PLAY;
      } else if (pp.state = PSTATE::LOOP_SILENCE) {
	cont.new_state = PSTATE::LOOP;
      }
      break;    
    case ' ':
      if(pp.state == PSTATE::PLAY){
	cont.cmd = COMMAND::LOOP_INIT;
	std::cout << "loop from: " << pp.offset << std::endl;
      } else if (pp.state == PSTATE::LOOP_REC) {
	cont.cmd = COMMAND::LOOP_FINISH;
	std::cout << "loop to: " << pp.offset << std::endl;
      }
      break;
    case 'b':
      cont.cmd = COMMAND::STATE_CHANGE;
      if(pp.state == PSTATE::PLAY){
	cont.new_state = PSTATE::BLOCK;
      } else if (pp.state == PSTATE::LOOP){
	cont.new_state = PSTATE::LOOP_BLOCK;
      } else if (pp.state == PSTATE::LOOP_BLOCK) {
	cont.new_state = PSTATE::LOOP;
	cv.notify_all();
      } else if (pp.state == PSTATE::BLOCK){
	cont.new_state = PSTATE::PLAY;
	cv.notify_all();
      }
      break;
    default:
      std::cout << "COMMAND NOT ACCEPTED!" << std::endl;
    }
    cmd_queue.push(cont);
  }

  // on exit, close stream  
  try {
    // Stop the stream
    dac.stopStream();
  } catch (RtAudioError &e) {
    e.printMessage();
  }
  
  if (dac.isStreamOpen()) {
    dac.closeStream();
  }
  
  return exit_code;
}
