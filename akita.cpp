#include <iostream>
#include <algorithm>
#include <iterator>
#include <sndfile.hh>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <boost/interprocess/ipc/message_queue.hpp>
#include "RtAudio.h"

namespace po = boost::program_options;
namespace lfree = boost::lockfree;
namespace iproc = boost::interprocess;

#define CHUNKSIZE 1024

// format enum
enum RWTYPES { UCHAR, SHORT, FLOAT, DOUBLE };

namespace PSTATE{
  enum PSTATE  { PLAY, LOOP_REC, LOOP, SILENCE };
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


struct play_params {
  // loop params
  lfree::spsc_queue<PSTATE::PSTATE>* cmd_queue;
  PSTATE::PSTATE state;
  bool loop;
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

// the parameterized callback function ...
template <typename READ_TYPE, typename WRITE_TYPE>
int abusive_play(void *outputBuffer, void *inputBuffer,
                 unsigned int nBufferFrames, double streamTime,
                 RtAudioStreamStatus status, void *userData) {

  // get the parameter container from the user data ...
  file_container<READ_TYPE> *fc = reinterpret_cast<file_container<READ_TYPE> *>(userData);
  params_to_abuse *params = fc->params;
  play_params *pp = fc->pp;
  lfree::spsc_queue<PSTATE::PSTATE> *cmd_queue = pp->cmd_queue;
  WRITE_TYPE *out_buf = (WRITE_TYPE *)outputBuffer;

  //get new play state ...
  PSTATE::PSTATE new_state;
  while(cmd_queue->pop(&new_state)){
      pp->state = new_state;
  }
  
  if (status) { std::cout << "Stream underflow detected!" << std::endl; }

  // SILENCE ! I KILL YOU !!
  if (pp->state == PSTATE::SILENCE) {
    for (int i = 0; i < nBufferFrames * params->buffer_cut; i++) {
      *out_buf++ = 0;
    }
    return 0;
  }
  
  // transfer samples from file buffer to output !
  for (int i = 0; i < nBufferFrames * params->buffer_cut; i++) {
    for (float j = 0; j < params->sample_repeat; j += 1.0) {
      if (pp->offset + i > fc->samples) {
        pp->offset = 0;
      }
      *out_buf++ = fc->file_buffer[pp->offset + i];
    }
  }

  // increament read offest

  pp->offset += nBufferFrames * params->offset_cut;

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
      ("cmd-queue", po::value<std::string>(), "Command Queue")
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

  bool use_stdout = true;
  std::string cmd_queue_name = "";
  iproc::message_queue* in_q;
  if (vm.count("cmd-queue")) {
    cmd_queue_name = vm["cmd-queue"].as<std::string>();
    if(cmd_queue_name.size() == 0){
      std::cout << "Please enter command queue name !" << std::endl;
      return 0;
    }
    use_stdout = false;
    in_q = new iproc::message_queue
      (iproc::open_or_create        //only create
       ,cmd_queue_name.c_str()       //name
       ,100
       ,sizeof(char)
       );

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

  lfree::spsc_queue<PSTATE::PSTATE> cmd_queue(10);
  
  // play params, for looping etc ...
  play_params pp;

  pp.state = PSTATE::PLAY; 
  pp.cmd_queue = &cmd_queue;
  pp.loop = false;
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
    bool get_input = true;
    while(get_input){
      if(use_stdout){
	std::cin.get(input);
      } else {
	unsigned int priority;
	iproc::message_queue::size_type recvd_size;
	in_q->receive(&input , sizeof(input ), recvd_size, priority);
      }

      switch(input) {
      case 'q':
	try {
	  // Stop the stream
	  dac.stopStream();
	} catch (RtAudioError &e) {
	  e.printMessage();
	}
	if (dac.isStreamOpen())
	  dac.closeStream();
	get_input = false;
	break;
      case 's':
	//PSTATE::PSTATE new_state = PSTATE::SILENCE;
	cmd_queue.push(PSTATE::SILENCE);
	break;
      case 'p':
	//PSTATE::PSTATE new_state = PSTATE::SILENCE;
	cmd_queue.push(PSTATE::PLAY);
	break;
      case 'o':
	std::cout << pp.offset << std::endl;
	break;
      default:
	std::cout << "COMMAND NOT ACCEPTED!" << std::endl;
      }
    }
  
  
  return exit_code;
}
