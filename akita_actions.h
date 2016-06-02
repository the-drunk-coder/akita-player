#include <condition_variable>
#include "akita_structures.h"
#include <iostream>

namespace akita_actions {

template <typename READ_TYPE>
void change_gain (source_params<READ_TYPE>& spar, filter_params& fpar, float gain_mod) {
  filter_command_container fcont(COMMAND::GAIN_CHANGE);
  fcont.gain = fpar.gain + gain_mod >= 1.0 ? 1.0 : fpar.gain + gain_mod;  
  fpar.cmd_queue->push(fcont);
 }

template <typename READ_TYPE>
void change_flippiness (source_params<READ_TYPE>& spar, filter_params& fpar, float flippiness_mod) {
  filter_command_container fcont(COMMAND::FLIPPINESS_CHANGE);
  fcont.flippiness = fpar.flippiness + flippiness_mod >= 1.0 ? 1.0 : fpar.flippiness + flippiness_mod;  
  fpar.cmd_queue->push(fcont);
 }

template <typename READ_TYPE>
void change_fuzziness (source_params<READ_TYPE>& spar, filter_params& fpar, float fuzziness_mod) {
  source_command_container scont(COMMAND::FUZZINESS_CHANGE);
  scont.new_fuzziness = spar.fuzziness + fuzziness_mod >= 1.0 ? 1.0 : spar.fuzziness + fuzziness_mod;  
  spar.cmd_queue->push(scont);
 }

template <typename READ_TYPE>
void change_samplerate (source_params<READ_TYPE>& spar, filter_params& fpar, int samplerate_mod) {
  source_command_container scont(COMMAND::SAMPLERATE_CHANGE);  
  scont.new_sample_repeat = spar.sample_repeat + samplerate_mod <= 0 ? 1 : spar.sample_repeat + samplerate_mod;  
  spar.cmd_queue->push(scont);
 }
 
template <typename READ_TYPE>
void toggle_filter (source_params<READ_TYPE>& spar, filter_params& fpar) {
  filter_command_container fcont(COMMAND::TOGGLE_FILTER);
  fpar.cmd_queue->push(fcont);
 }

template <typename READ_TYPE>
void toggle_reverb (source_params<READ_TYPE>& spar, filter_params& fpar) {
  filter_command_container fcont(COMMAND::TOGGLE_REVERB);
  fpar.cmd_queue->push(fcont);
 }

template <typename READ_TYPE>
void toggle_mean_filter (source_params<READ_TYPE>& spar, filter_params& fpar) {
  filter_command_container fcont(COMMAND::TOGGLE_MEAN_FILTER);
  fpar.cmd_queue->push(fcont);
 }

template <typename READ_TYPE>
void toggle_mute (source_params<READ_TYPE>& spar, filter_params& fpar) {
  source_command_container scont(COMMAND::STATE_CHANGE);
  if(spar.state == PLAY){
    scont.new_state = SILENCE;
  } else if (spar.state == LOOP){
    scont.new_state = LOOP_SILENCE;
  } else if (spar.state == SILENCE){
    scont.new_state = PLAY;
  } else if (spar.state == LOOP_SILENCE) {
    scont.new_state = LOOP;
  }
  spar.cmd_queue->push(scont);
 }
 
template <typename READ_TYPE>
void toggle_loop_state (source_params<READ_TYPE>& spar, filter_params& fpar) {  
  if(spar.state == PLAY){
    source_command_container scont(COMMAND::LOOP_INIT);
    std::cout << "loop from: " << spar.offset << " (~" << (float) spar.offset / spar.fc->samples  << ")" << std::endl;
    spar.cmd_queue->push(scont);
  } else if (spar.state == LOOP_REC) {
    source_command_container scont(COMMAND::LOOP_FINISH);
    std::cout << "loop to: " << spar.offset << " (~" << (float) spar.offset / spar.fc->samples  << ")" << std::endl;
    spar.cmd_queue->push(scont);
  } else if (spar.state == LOOP) {
    source_command_container scont(COMMAND::LOOP_RELEASE);
    std::cout << "loop off" << std::endl;
    spar.cmd_queue->push(scont);
  }
 }

template <typename READ_TYPE>
void toggle_block (source_params<READ_TYPE>& spar, filter_params& fpar, std::condition_variable& cv) {  
  source_command_container scont(COMMAND::STATE_CHANGE);
  if(spar.state == PLAY){
    scont.new_state = BLOCK;
  } else if (spar.state == LOOP){    
    scont.new_state = LOOP_BLOCK;
  } else if (spar.state == LOOP_BLOCK) {    
    scont.new_state = LOOP;
    cv.notify_all();
  } else if (spar.state == BLOCK){    
    scont.new_state = PLAY;
    cv.notify_all();
  }
  spar.cmd_queue->push(scont);
}

template <typename READ_TYPE>
void toggle_filter_band (source_params<READ_TYPE>& spar, filter_params& fpar, char input){
  // stupidly simple, but it seems to work ...
  if (input == '0'){input = 9;}
  else {input -= 49;}
  fpar.fbank->toggle_band(input);  
 }
 
}
