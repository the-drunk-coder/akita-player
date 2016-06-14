#ifndef _AKITA_STRUCTURES_H_
#define _AKITA_STRUCTURES_H_

#include <iostream>
#include <map>
#include <boost/algorithm/string.hpp>
#include "RtAudio.h"
#include <sndfile.hh>
#include "akita_filters.h"
#include <stk/FreeVerb.h>
#include <memory>
#include <atomic>

// format enum
enum INTERFACE {PLAIN, ADVANCED, OSC};

// custom stream to extract enum from command line
std::istream &operator>>(std::istream &in, INTERFACE &iface) {
  std::string token;
  in >> token;

  boost::to_upper(token);

  if (token == "PLAIN") {
    iface = PLAIN;
  } else if (token == "ADVANCED") {
    iface = ADVANCED;
  } else if (token == "OSC") {
    iface = OSC;
  } 
  return in;
}

// format enum
enum RWTYPES { UCHAR, SHORT, FLOAT, DOUBLE };

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
enum PSTATE { WAIT_EVENT_START, PLAY, PLAY_EVENT, LOOP_REC, LOOP, LOOP_SILENCE, LOOP_BLOCK, SILENCE, BLOCK };

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

struct akita_play_event {
  enum STATE { NEW, IN_PROGRESS, FINISHED };

  STATE state = NEW;
    
  float start;
  int dur;
  int length_samples;
  int samples_played = 0;
  
  akita_play_event () {
    start = 0.0;
    dur = 0;
    state = FINISHED;
  }
  
  akita_play_event (float start, int dur) {
    this->start = start;
    this->dur = dur;
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

  void update_range (akita_play_event& ev) {
    start_sample = (frames * ev.start) * channels;
    end_sample = start_sample + (((samplerate / 1000) * ev.dur) * channels);
    if(end_sample >= frames){
      end_sample = frames;
    }
    ev.length_samples = end_sample - start_sample;
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
 * container for command line options
 */
struct options_container {

  INTERFACE iface;
  
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
  
  //int channel_offset;
  int out_channels;
  float pan;

  int mean_filter_points;

  float start;
  float end;

  float reverb_mix;
  
  bool mono = true;

  int udp_port;
};

/*
 * Parameter structs to pass to noise- and filter threads.
 */
// params to pass to noise source
template <typename READ_TYPE>
struct source_params {
  
  file_container<READ_TYPE> fc;

  // current event, only necessary for osc mode 
  std::unique_ptr<akita_play_event> current_event;

  INTERFACE iface;
  
  // state ... playing, silent, blocking, looping etc
  std::atomic<PSTATE> state;

  // loop params
  std::atomic<long int> loop_start;
  std::atomic<long int> loop_end;

  // the current position within the file buffer
  std::atomic<long int> offset;

  // buffer glitch parameters
  std::atomic<float> buffer_cut;
  std::atomic<float> offset_cut;
  std::atomic<float> sample_repeat;

  std::atomic<float> fuzziness;
  
  source_params(options_container& opts) :     
     fc(opts.filename, opts.start, opts.end, opts.mono),     
     current_event(new akita_play_event())  
  {
    iface = opts.iface;
    if(opts.iface == OSC){
      state = SILENCE;
    } else {
      state = PLAY;
    }
    loop_start = 0;
    loop_start = 0;
    offset = 0;
    
    state = opts.initial_state;
    sample_repeat = opts.sample_repeat;
    fuzziness = opts.fuzziness;
        
    // set starting point
    offset = fc.start_sample;
    
    // default handling ...
    if(opts.offset_cut == 2){      
      offset_cut = fc.channels;
      opts.offset_cut = offset_cut;
    } else {
      offset_cut = opts.offset_cut;
    }

    if(opts.buffer_cut == 2){      
      buffer_cut = fc.channels;
      opts.buffer_cut = buffer_cut;
    } else {
      buffer_cut = opts.buffer_cut;
    }
               
  }  
};

// params for the filter 
struct filter_params {
  
  filterbank fbank; 
  mean_filterbank m_fbank;
  canonical_sos_filter lowpass;
  stk::FreeVerb rev;
  
  PMODE mode;
  int channels;
  
  std::atomic<float> gain;
  
  std::atomic<bool> lowpass_on;
  std::atomic<bool> filterbank_on;
  std::atomic<bool> mean_filter_on;
  std::atomic<bool> reverb_on;
  
  std::atomic<float> flippiness;

  // pan
  std::atomic<int> pan_offset;
  std::atomic<float> pan_ratio;

  float* frame_buffer;

  void update_pan(float pan){    
    pan_offset = (int) pan;
    pan_ratio = pan - pan_offset;
  }
  
  filter_params (options_container& opts) :
    fbank(opts.out_channels, opts.samplerate, 100, 8000, 10),
    m_fbank(opts.out_channels, 5),
    lowpass()
  {

    lowpass_on = false;
    filterbank_on = false;
    reverb_on = false;
    mean_filter_on = false;
    
    mode = opts.initial_mode;
    gain = opts.initial_gain;

    update_pan(opts.pan);
    
    flippiness = opts.flip_prob;
    this->channels = opts.out_channels;
    frame_buffer = new float[channels];
                
    if(opts.mean_filter_points > 0){
      m_fbank.update(opts.mean_filter_points);
      mean_filter_on = true;
    }     
  }

  ~filter_params () {
    //std::cout << "cleaning filter params" << std::endl;    
    delete frame_buffer;
  }
};

#endif /* _AKITA_STRUCTURES_H_*/
