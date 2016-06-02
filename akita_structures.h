#ifndef _AKITA_STRUCTURES_H_
#define _AKITA_STRUCTURES_H_

#include <iostream>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include "RtAudio.h"
#include <sndfile.hh>
#include "akita_filters.h"
#include <stk/FreeVerb.h>

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
  enum COMMAND { MODE_CHANGE,
		 STATE_CHANGE,
		 SAMPLERATE_CHANGE,
		 FUZZINESS_CHANGE,
		 FLIPPINESS_CHANGE,
		 GAIN_CHANGE,
		 TOGGLE_MEAN_FILTER,
		 TOGGLE_FILTER,
		 TOGGLE_REVERB,
		 LOOP_INIT,
		 LOOP_FINISH,
		 LOOP_RELEASE
  };
}

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

  source_command_container(){};
  
  source_command_container(COMMAND::COMMAND cmd){
    this->cmd = cmd;
  }
};

struct filter_command_container {
  COMMAND::COMMAND cmd;
  
  // parameters that could be subject to change
  PMODE new_mode;
  float gain;
  float flippiness;

  filter_command_container(){};
  
  filter_command_container(COMMAND::COMMAND cmd){
    this->cmd = cmd;
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
  bool reverb = false;
  
  float flippiness;
  
  filter_params (options_container& opts, int channels) {
    mode = opts.initial_mode;
    gain = opts.initial_gain;
    flippiness = opts.flip_prob;
    this->channels = channels;
    
    cmd_queue = new lfree::spsc_queue<filter_command_container>(10);
    fbank = new filterbank(channels, opts.samplerate, 100, 8000, 10);
    
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
  }
};

#endif /* _AKITA_STRUCTURES_H_*/
