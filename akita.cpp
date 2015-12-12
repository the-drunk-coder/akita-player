#include <iostream>
#include <algorithm>
#include <iterator>
#include <sndfile.hh>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include "RtAudio.h"

namespace po = boost::program_options;

#define CHUNKSIZE 1024

// format enum
enum RWTYPES {
  INT8,
  INT16,
  FLOAT32,
  FLOAT64
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

// file container to communicate with callback function
template<typename READ_TYPE>
class file_container {
public:
  file_container(SndfileHandle file){
    samples = file.frames() * file.channels();
    file_buffer = new READ_TYPE[samples];

    offset = 0;
    channels = file.channels();
  }

  ~file_container(){
    delete[] file_buffer;
  }

  READ_TYPE* file_buffer;
  long int offset;
  long int samples;
  short channels;

  // the params to make things weird!
  params_to_abuse* params;
};


// the parameterized callback function ...
template<typename READ_TYPE, typename WRITE_TYPE>
int abusive_play( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
         double streamTime, RtAudioStreamStatus status, void *userData )
{
  // get the parameter container from the user data ...
  file_container<READ_TYPE>* fc = reinterpret_cast<file_container<READ_TYPE>*>(userData);
  params_to_abuse* params = fc->params;
  WRITE_TYPE* out_buf = (WRITE_TYPE*) outputBuffer;

  if ( status ) { std::cout << "Stream underflow detected!" << std::endl; }

  // transfer samples from file buffer to output !
  for(int i = 0; i < nBufferFrames * params->buffer_cut; i++) {
    for(float j = 0; j < params->sample_repeat; j+=1.0) {
      if(fc->offset + i > fc->samples) {
        fc->offset = 0;
      }
      *out_buf++ = fc->file_buffer[fc->offset + i];
    }
  }

  // increament read offest
  fc->offset += nBufferFrames * params->offset_cut;
  return 0;
}

/*
* INITIALIZE COMMAND LINE OPTIONS !!
*/

// A helper function to simplify the main part (from po example ...)
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}

// custom validator for r/w types
std::istream& operator>> (std::istream &in, RWTYPES& rwtype)
{
    std::string token;
    in >> token;

    boost::to_upper(token);

    if (token.compare("INT8")) {
        rwtype = RWTYPES::INT8;
    } else if (token.compare("INT16")) {
      rwtype = RWTYPES::INT16;
    } else if (token.compare("FLOAT32")) {
        rwtype = RWTYPES::FLOAT32;
    } else if (token.compare("FLOAT64")) {
        rwtype = RWTYPES::FLOAT64;
    }

    return in;
}

// initialize the command line parameter parser
po::options_description init_opts(int ac, char* av[], po::variables_map* vm, params_to_abuse *params){
  po::options_description desc("Parameters to use in a creative way");
  desc.add_options()
      ("help", "Display this help!")
      ("input-file", po::value<std::string>(),
            "The input file - WAV!")
      ("sample-repeat", po::value<float>(),
            "Repeat every sample n times!")
      ("buffer-cut", po::value<float>(),
            "Don't fill output buffer completely!")
      ("offset-cut", po::value<float>(),
            "Modify offset increment (chunk size read from buffer)!")
      ("read-type", po::value<RWTYPES>(&(*params).read_type),
            "Type used to read from audio file!")
      ("write-type", po::value<RWTYPES>(&(*params).write_type),
            "Type used to write to audio buffer!")
      ("stream-type", po::value<RWTYPES>(&(*params).stream_type),
            "Type used for the audio stream!")
  ;
  // ----- end options ... what kind of syntax is this ??

  po::positional_options_description p;
  p.add("input-file", -1);

  po::store(po::command_line_parser(ac, av).
            options(desc).positional(p).run(), *vm);
  po::notify(*vm);
  return desc;
}

template<typename READ_TYPE, typename WRITE_TYPE>
int load_file_and_play(SndfileHandle file, params_to_abuse *params) {

  file_container<READ_TYPE> *fc = new file_container<READ_TYPE>(file);
  fc->params = params;

  int frame_chunksize = CHUNKSIZE * file.channels();
  READ_TYPE *chunk_buffer = new READ_TYPE[frame_chunksize];

  // Read file into buffer
  long int read_offset = 0;
  while(file.readf(chunk_buffer, CHUNKSIZE) == CHUNKSIZE){
    for (int i = 0; i < frame_chunksize; i++){
      fc->file_buffer[read_offset + i] = chunk_buffer[i];
    }
    read_offset += frame_chunksize;
  }
  // not needed any longer !
  delete[] chunk_buffer;

	RtAudio dac;
  if ( dac.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    return 0;
  }
  RtAudio::StreamParameters parameters;
  parameters.deviceId = dac.getDefaultOutputDevice();
  parameters.nChannels = 2;
  parameters.firstChannel = 0;
  unsigned int sampleRate = 44100;
  unsigned int bufferFrames = 1024; // 256 sample frames

  try {
    if(params->stream_type == RWTYPES::INT8) {
      dac.openStream( &parameters, NULL, RTAUDIO_SINT8,
                      sampleRate, &bufferFrames, &abusive_play<READ_TYPE, WRITE_TYPE>, (void*) fc);
    } else if (params->stream_type == RWTYPES::FLOAT32) {
      dac.openStream( &parameters, NULL, RTAUDIO_FLOAT32,
                      sampleRate, &bufferFrames, &abusive_play<READ_TYPE, WRITE_TYPE>, (void*) fc);
    } else if (params->stream_type == RWTYPES::FLOAT64) {
      dac.openStream( &parameters, NULL, RTAUDIO_FLOAT64,
                      sampleRate, &bufferFrames, &abusive_play<READ_TYPE, WRITE_TYPE>, (void*) fc);
    } else {
      // this is probably the most commom option
      dac.openStream( &parameters, NULL, RTAUDIO_SINT16,
                      sampleRate, &bufferFrames, &abusive_play<READ_TYPE, WRITE_TYPE>, (void*) fc);
    }
    dac.startStream();
  } catch ( RtAudioError& e ) {
    e.printMessage();
    return 0;
  }

  char input;
  std::cout << "\nPlaying ... press <enter> to quit.\n";
  std::cin.get( input );
  try {
    // Stop the stream
    dac.stopStream();
  } catch (RtAudioError& e) {
    e.printMessage();
  }
  if ( dac.isStreamOpen() ) dac.closeStream();

  return 0;
}

/*
* MAIN ROUTINE!!
*/
int main (int ac, char* av[]) {
  // Salutations!
  std::cout << "\n~~ akita - create noise abusing low-level audio parameters! ~~\n"<< std::endl ;

  // get params to abuse ...
  params_to_abuse params;

  // initialize command line options
  po::variables_map vm;
  po::options_description desc = init_opts(ac, av, &vm, &params);

  // help output
  if (vm.count("help")) {
      std::cout << "usage: akita <file> [options] \n"<< std::endl ;
      std::cout << desc;
      return 0;
  }

  // get filename to use ...
  const char * fname = "" ;
  if (vm.count("input-file")) {
        fname = vm["input-file"].as<std::string>().c_str();
  } else {
    std::cout << "Please specify input file !" << std::endl;
    return 0;
  }

  SndfileHandle file = SndfileHandle (fname) ;

  params.sample_repeat = 1.0;
  if (vm.count("sample-repeat")) {
        params.sample_repeat = vm["sample-repeat"].as<float>();
  }

  params.buffer_cut = (float) file.channels();
  if (vm.count("buffer-cut")) {
        params.buffer_cut = vm["buffer-cut"].as<float>();
  }

  params.offset_cut = (float) file.channels();
  if (vm.count("offset-cut")) {
        params.offset_cut = vm["offset-cut"].as<float>();
  }

  params.stream_type = RWTYPES::INT16;
  if (vm.count("stream-type")) {
        params.stream_type = vm["stream-type"].as<RWTYPES>();
  }

  params.read_type = RWTYPES::INT16;
  if (vm.count("stream-type")) {
        params.read_type = vm["read-type"].as<RWTYPES>();
  }

  params.write_type = RWTYPES::INT16;
  if (vm.count("stream-type")) {
        params.read_type = vm["write-type"].as<RWTYPES>();
  }

  // display some file info ...
  std::cout << "Input file: " << fname << ", " << file.channels() << "ch, "
    << file.samplerate() << "Hz, "<< file.frames() << " frames." << std::endl;

  // ok, here it get's a little awkward ... a dynamically-typed language would come
  // in handy here ...
  int exit_code = 0;

  if (params.read_type == RWTYPES::FLOAT32) {
    if(params.write_type == RWTYPES::INT8) {
      exit_code = load_file_and_play<float_t, int8_t>(file, &params);
    } else if(params.write_type == RWTYPES::FLOAT32) {
      exit_code = load_file_and_play<float_t, float_t>(file, &params);
    } else if(params.write_type == RWTYPES::FLOAT64) {
      exit_code = load_file_and_play<float_t, double_t>(file, &params);
    } else {
      exit_code = load_file_and_play<float_t, int16_t>(file, &params);
    }
  } else if (params.read_type == RWTYPES::FLOAT64) {
    if(params.write_type == RWTYPES::INT8) {
      exit_code = load_file_and_play<double_t, int8_t>(file, &params);
    } else if(params.write_type == RWTYPES::FLOAT32) {
      exit_code = load_file_and_play<double_t, float_t>(file, &params);
    } else if(params.write_type == RWTYPES::FLOAT64) {
      exit_code = load_file_and_play<double_t, double_t>(file, &params);
    } else {
      exit_code = load_file_and_play<double_t, int16_t>(file, &params);
    }
  } else {
    if(params.write_type == RWTYPES::INT8) {
      exit_code = load_file_and_play<int16_t, int8_t>(file, &params);
    } else if(params.write_type == RWTYPES::FLOAT32) {
      exit_code = load_file_and_play<int16_t, float_t>(file, &params);
    } else if(params.write_type == RWTYPES::FLOAT64) {
      exit_code = load_file_and_play<int16_t, double_t>(file, &params);
    } else {
      exit_code = load_file_and_play<int16_t, int16_t>(file, &params);
    }
  }

  return exit_code;
}
