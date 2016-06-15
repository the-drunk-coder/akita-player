#ifndef _AKITA_ACTIONS_H_
#define _AKITA_ACTIONS_H_

#include <condition_variable>
#include "akita_structures.h"
#include <iostream>

namespace akita_actions {

template <typename READ_TYPE>
void change_gain (source_params<READ_TYPE>& spar, filter_params& fpar, float new_gain) {  
  if (new_gain < 0) { fpar.gain = 0; }
  else if (new_gain > 1.0) { fpar.gain = 1.0; }
  else { fpar.gain = new_gain; }
  std::cout << "gain change, new gain: " << fpar.gain << std::endl;
 }

template <typename READ_TYPE>
void change_reverb_mix (source_params<READ_TYPE>& spar, filter_params& fpar, float new_reverb_mix) {  
  fpar.rev.setEffectMix(new_reverb_mix);
  fpar.reverb_on = new_reverb_mix > 0.0;
 }

template <typename READ_TYPE>
void change_lowpass (source_params<READ_TYPE>& spar, filter_params& fpar, float new_q, float new_freq) {  
  fpar.lowpass.q = new_q;
  fpar.lowpass.frequency = new_freq;
  fpar.lowpass.update();
  fpar.lowpass_on = new_freq < 20000;
}

template <typename READ_TYPE>
void change_pan (source_params<READ_TYPE>& spar, filter_params& fpar, float new_pan) {  
  if(new_pan < 0) { fpar.update_pan(fpar.channels); }
  if(new_pan > fpar.channels) { fpar.update_pan(0.0); }
  else { fpar.update_pan(new_pan); }
}
 
template <typename READ_TYPE>
void change_flippiness (source_params<READ_TYPE>& spar, filter_params& fpar, float new_flippiness) {
  if (new_flippiness < 0) { fpar.flippiness = 0; }
  else if (new_flippiness > 1.0) { fpar.flippiness = 1.0; }
  else { fpar.flippiness = new_flippiness; }
  std::cout << "flippiness change, new flippiness: " << fpar.flippiness << std::endl;  
 }

template <typename READ_TYPE>
void change_fuzziness (source_params<READ_TYPE>& spar, filter_params& fpar, float new_fuzziness) {
  if (new_fuzziness < 0) { spar.fuzziness = 0; }
  else if (new_fuzziness > 1.0) { spar.fuzziness = 1.0; }
  else { spar.fuzziness = new_fuzziness; }
  std::cout << "fuzziness change, new fuzziness: " << spar.fuzziness << std::endl;    
 }

template <typename READ_TYPE>
void change_samplerate (source_params<READ_TYPE>& spar, filter_params& fpar, int new_samplerate) {
  if (new_samplerate < 1) { spar.sample_repeat = 1; }  
  else { spar.sample_repeat = new_samplerate; }
  std::cout << "samplerate change, new samplerate: " << spar.sample_repeat << std::endl;    
 }
 
template <typename READ_TYPE>
void toggle_filterbank (source_params<READ_TYPE>& spar, filter_params& fpar) {
  fpar.filterbank_on = ! fpar.filterbank_on;
  std::cout << "filter toggled, now: " << !fpar.filterbank_on << std::endl;
 }

template <typename READ_TYPE>
void toggle_reverb (source_params<READ_TYPE>& spar, filter_params& fpar) {
  fpar.reverb_on = !fpar.reverb_on;
  std::cout << "reverb toggled, now: " << !fpar.reverb_on << std::endl;
 }

template <typename READ_TYPE>
void toggle_mean_filter (source_params<READ_TYPE>& spar, filter_params& fpar) {
  fpar.mean_filter_on = !fpar.mean_filter_on;
  std::cout << "mean filter toggled, now: " << !fpar.mean_filter_on << std::endl;  
 }

template <typename READ_TYPE>
void toggle_lowpass_filter (source_params<READ_TYPE>& spar, filter_params& fpar) {
  fpar.lowpass_on = !fpar.lowpass_on;
  std::cout << "lowpass filter toggled, now: " << !fpar.mean_filter_on << std::endl;  
 }

template <typename READ_TYPE>
void toggle_mute (source_params<READ_TYPE>& spar, filter_params& fpar) {
  if(spar.state == PLAY){
    spar.state = SILENCE;
  } else if (spar.state == LOOP){
    spar.state= LOOP_SILENCE;
  } else if (spar.state == SILENCE){
    spar.state = PLAY;
  } else if (spar.state == LOOP_SILENCE) {
    spar.state = LOOP;
  }
  std::cout << "mute toggled, now: " << spar.state << std::endl;
 }
 
template <typename READ_TYPE>
void toggle_loop_state (source_params<READ_TYPE>& spar, filter_params& fpar) {  
  if(spar.state == PLAY){
    spar.state = LOOP_REC;
    std::cout << "loop from: " << spar.offset << " (~" << (float) spar.offset / spar.fc.samples  << ")" << std::endl;
  } else if (spar.state == LOOP_REC) {
    spar.state = LOOP;
    std::cout << "loop to: " << spar.offset << " (~" << (float) spar.offset / spar.fc.samples  << ")" << std::endl;    
  } else if (spar.state == LOOP) {
    spar.state = PLAY;
    std::cout << "loop off" << std::endl;    
  }
 }

template <typename READ_TYPE>
void toggle_block (source_params<READ_TYPE>& spar, filter_params& fpar, std::condition_variable& cv) {  
  if(spar.state == PLAY){
    spar.state = BLOCK;
  } else if (spar.state == LOOP){    
    spar.state = LOOP_BLOCK;
  } else if (spar.state == LOOP_BLOCK) {    
    spar.state = LOOP;
    cv.notify_all();
  } else if (spar.state == BLOCK){    
    spar.state = PLAY;
    cv.notify_all();
  }
  std::cout << "blocking toggled, now: " << spar.state << std::endl;  
}

template <typename READ_TYPE>
void toggle_filter_band (source_params<READ_TYPE>& spar, filter_params& fpar, char input){
  // stupidly simple, but it seems to work ...
  if (input == '0'){input = 9;}
  else {input -= 49;}
  fpar.fbank.toggle_band(input);
  fpar.fbank.print_bands();
 }
 
}

#endif
