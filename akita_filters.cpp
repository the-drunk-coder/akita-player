#include <cmath>
#include <iostream>
#include <cstdlib>
#include <climits>
#include <bitset>
#include "akita_filters.h"

float random_flip(float sample, float prob){

  float_conv sample_conv;
  sample_conv.float_val = sample;

  //use a bitset
  std::bitset<32> sample_bits(sample_conv.uint_val);

  if(prob >= 1.0) {
    sample_bits.flip();
  } else {
    for(int i = 0; i < 32; i++){
      if((((float) rand()) / INT_MAX) < prob){
	sample_bits.flip(i); 
      }
    }
  }
  
  // convert back
  sample_conv.uint_val = sample_bits.to_ulong();
  
  return sample_conv.float_val;
}

//canonicalsos filter methods
canonical_sos_filter::canonical_sos_filter () {
  update(5000, 2, 44100, LP);
}

canonical_sos_filter::canonical_sos_filter (float frequency, float q, int samplerate, FMODE mode) {
  update(frequency, q, samplerate, mode);
}

void canonical_sos_filter::update (float frequency, float q, int samplerate, FMODE mode) {    
  del1 = 0;
  del2 = 0;
  
  k = tanh( (M_PI * frequency) / samplerate);

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

float canonical_sos_filter::calculate (float sample) {
  float intermediate = sample + ((-1.0 * a1) * del1) + ((-1.0 * a2) * del2);
  float out = (b0 * intermediate) + (b1 * del1) + (b2 * del2);
  del2 = del1;
  del1 = intermediate;
  return out;
}

void canonical_sos_filter::process(float& sample){
  sample = calculate(sample);          
}

simple_mean_filter::simple_mean_filter(int points){
  init(points);
}

// initialize with 7 points per default ...
simple_mean_filter::simple_mean_filter () {
  init(13);
}
  
simple_mean_filter::~simple_mean_filter () {
  delete[] delay;
}

void simple_mean_filter::init (int points){
  if(initialized){
    delete[] delay;
  }
  
  this->points = points;
  newest = 0;
  //kernel = new float[points];
  factor = 1.0f / points;
  
  delay = new float[points];
  
  for (int i = 0; i < points; i++) {
    delay[i] = 0.0f;
  }
  
  initialized = true;
}

float simple_mean_filter::calculate (float sample) {
  float out = 0.0;
  
  newest++;
  
  if (newest >= points) {
    newest = 0;
  }
  
  delay[newest] = sample;
  
  // calculate convolution
  short delay_index = newest;
  for (int i = 0; i < points; i++) {
    
    //out += kernel[i] * delay[delay_index];
    out += factor * delay[delay_index];
    
    delay_index--;
    
    if (delay_index <= 0) {
      delay_index = points - 1;
    }
  }
  
  return out;        
}

void simple_mean_filter::apply(float& sample){
  sample = calculate(sample);
}

mean_filterbank::mean_filterbank (int channels, int points) {
  filters = new simple_mean_filter[channels];
  for(int i = 0; i < channels; i++){
    filters[i].init(points);
  }
}

mean_filterbank::~mean_filterbank () {
  delete[] filters;
}

void mean_filterbank::apply (int channel, float& sample) {
  filters[channel].apply(sample);
}

filterbank::filterbank(int channels, int samplerate, float lowcut, float hicut, int bands) {
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

filterbank::~filterbank () {
  delete [] fbank;
  delete [] fmask;
}

void filterbank::apply (int channel, float& sample) {    
  for(int b = 0; b < bands; b++){
    if(fmask[b]){	
      fbank[(channel * bands) + b].process(sample);
    }
  }
}

void filterbank::toggle_band (int band) {
  if(fmask[band]) {
    fmask[band] = false;
  } else {
    fmask[band] = true;
  } 
}

void filterbank::print_bands () {
  std::cout << "Filter bands: ";
  for(int i = 0; i < bands; i++){
    std::cout << "[" << fmask[i] << "] ";
  }
  std::cout << std::endl;
}