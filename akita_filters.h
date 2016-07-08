#ifndef _AKITA_FILTERS_H_
#define _AKITA_FILTERS_H_

#include <cstdint>

/*
 * random bit flipping
 */
union float_conv {
  uint32_t uint_val;
  float float_val;
};

float random_flip(float sample, float prob);

struct peak_filter {

  float frequency;
  float bandwidth;
  float gain;
  float samplerate;
  
  // parameters needed for peak filtering 
  float h_zero;
  float v_zero;
  float d;

  // delay samples
  float del1, del2;

  // boost/cut factor
  float c;
  
  peak_filter();

  void update();
  void update(float frequency, float gain, float bandwith);

  float calculate(float sample);
  void process(float& sample);
};

/*
 * Filter and filterbank
 */
struct canonical_sos_filter {
  enum FMODE {HP, BP, LP, NOTCH, AP};

  float frequency;
  float q;
  float gain;
  FMODE mode;
  int samplerate;
  
  float a1, a2;
  float b0, b1, b2;
  float k;  
  
  float del1, del2; 

  
  canonical_sos_filter();
  canonical_sos_filter (FMODE mode);
  canonical_sos_filter (float frequency, float q, float gain, int samplerate, FMODE mode);

  // update filter parameters
  void update ();
  void update (float frequency, float q, float gain, int samplerate, FMODE mode);
  
  // calculate filtered sample
  float calculate (float sample);

  // process sample,for pipeline usage
  void process (float& sample);
};


/*
 * simple mean filter, to shave the edge off a little ...
 */
struct simple_mean_filter {
  //float* kernel;
  float* delay;
  float factor;

  bool initialized = false;
  
  short newest;
  
  short points;

  simple_mean_filter (int points);

  // initialize with 13 points per default ...
  simple_mean_filter ();
  
  ~simple_mean_filter ();

  void init (int points);
  
  float calculate (float sample);
  
  void apply (float& sample);

};

struct mean_filterbank {
  int channels;
  simple_mean_filter* filters;
  
  mean_filterbank (int channels, int points);

  ~mean_filterbank();
  
  void update(int points);
  
  void apply (int channel, float& sample);
};


// a simple filterbank consisting of several sos filters
struct filterbank {
  
  int bands;
  int channels;
  canonical_sos_filter* fbank;
  bool* fmask;

  filterbank(int channels, int samplerate, float lowcut, float hicut, int bands);

  ~filterbank();
  
  void apply (int channel, float& sample);

  void toggle_band (int band);

  void print_bands ();
};

#endif
