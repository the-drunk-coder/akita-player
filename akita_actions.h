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
 }

template <typename READ_TYPE>
void change_reverb_mix (source_params<READ_TYPE>& spar, filter_params& fpar, float new_reverb_mix) {  
  fpar.rev.setEffectMix(new_reverb_mix);
  fpar.reverb_on = new_reverb_mix > 0.0;
 }

// first lti ops 
template <typename READ_TYPE>
void change_hipass (source_params<READ_TYPE>& spar, filter_params& fpar, float new_q, float new_freq) {  
  if(new_q < 0.1) { new_q = 0.1; }
  if(new_freq < 1) { new_freq = 1; }
  if(new_freq > 20000) { new_freq = 20000; }
  fpar.hipass.q = new_q;
  fpar.hipass.frequency = new_freq;
  fpar.hipass.update_internals();
  fpar.hipass_on = new_freq > 2;
}

template <typename READ_TYPE>
  void change_peak (source_params<READ_TYPE>& spar, filter_params& fpar, float new_bandwidth, float new_freq, float new_gain) {  
  if(new_bandwidth < 0) { new_bandwidth = 1; }
  if(new_freq < 1) { new_freq = 1; }
  if(new_freq > 20000) { new_freq = 20000; }
  //if(new_gain > 1.1) {new_gain = 1.1;}
  //if(new_gain < -1.0) {new_gain = -1.0;}
  fpar.peak.bandwidth = new_bandwidth;
  fpar.peak.frequency = new_freq;
  fpar.peak.gain = new_gain;
  fpar.peak.update_internals();
  fpar.peak_on = new_gain != 0.0 && new_gain != -0.0;  
}

template <typename READ_TYPE>
void change_lowpass (source_params<READ_TYPE>& spar, filter_params& fpar, float new_q, float new_freq) {  
  if(new_q < 0.1) { new_q = 0.1; }
  if(new_freq < 1) { new_freq = 1; }
  if(new_freq > 20000) { new_freq = 20000; }
  fpar.lowpass.q = new_q;
  fpar.lowpass.frequency = new_freq;
  fpar.lowpass.update_internals();
  fpar.lowpass_on = new_freq < 20000;
}

template <typename READ_TYPE>
  void change_nonlin (source_params<READ_TYPE>& spar, filter_params& fpar,
		      float kn, float kp, float gn,
		      float gp, float alpha_mix, float gain_sc,
		      float g_pre, float g_post, float lp_freq, bool on) {
  fpar.nonlin.kn = kn;
  fpar.nonlin.kp = kp;
  fpar.nonlin.gn = gn;
  fpar.nonlin.gp = gp;
  fpar.nonlin.alpha_mix = alpha_mix;
  fpar.nonlin.gain_sc = gain_sc;
  fpar.nonlin.g_pre = g_pre;
  fpar.nonlin.g_post = g_post;
  fpar.nonlin.lp_freq = lp_freq;
  fpar.nonlin.update_internals();
  fpar.nonlin_on = on;
}

// second lti ops 
template <typename READ_TYPE>
void change_hipass_2 (source_params<READ_TYPE>& spar, filter_params& fpar, float new_q, float new_freq) {  
  if(new_q < 0.1) { new_q = 0.1; }
  if(new_freq < 1) { new_freq = 1; }
  if(new_freq > 20000) { new_freq = 20000; }
  fpar.hipass_2.q = new_q;
  fpar.hipass_2.frequency = new_freq;
  fpar.hipass_2.update_internals();
  fpar.hipass_2_on = new_freq > 2;
}

template <typename READ_TYPE>
  void change_peak_2 (source_params<READ_TYPE>& spar, filter_params& fpar, float new_bandwidth, float new_freq, float new_gain) {  
  if(new_bandwidth < 0) { new_bandwidth = 1; }
  if(new_freq < 1) { new_freq = 1; }
  if(new_freq > 20000) { new_freq = 20000; }
  //if(new_gain > 1.1) {new_gain = 1.1;}
  //if(new_gain < -1.0) {new_gain = -1.0;}
  fpar.peak_2.bandwidth = new_bandwidth;
  fpar.peak_2.frequency = new_freq;
  fpar.peak_2.gain = new_gain;
  fpar.peak_2.update_internals();
  fpar.peak_2_on = new_gain != 0.0 && new_gain != -0.0;  
}

template <typename READ_TYPE>
void change_lowpass_2 (source_params<READ_TYPE>& spar, filter_params& fpar, float new_q, float new_freq) {  
  if(new_q < 0.1) { new_q = 0.1; }
  if(new_freq < 1) { new_freq = 1; }
  if(new_freq > 20000) { new_freq = 20000; }
  fpar.lowpass_2.q = new_q;
  fpar.lowpass_2.frequency = new_freq;
  fpar.lowpass_2.update_internals();
  fpar.lowpass_2_on = new_freq < 20000;
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
 }

template <typename READ_TYPE>
void change_fuzziness (source_params<READ_TYPE>& spar, filter_params& fpar, float new_fuzziness) {
  if (new_fuzziness < 0) { spar.fuzziness = 0; }
  else if (new_fuzziness > 1.0) { spar.fuzziness = 1.0; }
  else { spar.fuzziness = new_fuzziness; }
 }

template <typename READ_TYPE>
void change_sample_repeat (source_params<READ_TYPE>& spar, filter_params& fpar, int new_sample_repeat) {
  if (new_sample_repeat < 1) { spar.sample_repeat = 1; }  
  else { spar.sample_repeat = new_sample_repeat; }  
 }

template <typename READ_TYPE>
void change_samplerate_mod (source_params<READ_TYPE>& spar, filter_params& fpar, float new_samplerate_mod) {
  if (new_samplerate_mod <= 0.1) { spar.samplerate_mod = 0.1; }  
  else { spar.samplerate_mod = new_samplerate_mod; }
 }

template <typename READ_TYPE>
void toggle_filterbank (source_params<READ_TYPE>& spar, filter_params& fpar) {
  fpar.filterbank_on = ! fpar.filterbank_on;
 }

template <typename READ_TYPE>
void toggle_reverb (source_params<READ_TYPE>& spar, filter_params& fpar) {
  fpar.reverb_on = !fpar.reverb_on;
 }

template <typename READ_TYPE>
void toggle_mean_filter (source_params<READ_TYPE>& spar, filter_params& fpar) {
  fpar.mean_filter_on = !fpar.mean_filter_on;
 }

template <typename READ_TYPE>
void toggle_lowpass_filter (source_params<READ_TYPE>& spar, filter_params& fpar) {
  fpar.lowpass_on = !fpar.lowpass_on;
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
