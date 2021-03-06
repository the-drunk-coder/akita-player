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


struct filter {
  
  void process(float& sample);
  virtual float calculate(float sample) = 0;
  virtual void update_internals() = 0;
  
};

struct peak_filter : public filter {

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

  void update_internals();
  float calculate(float sample);  
};

/*
 * Filter and filterbank
 */
struct canonical_sos_filter : public filter {
  enum FMODE {HP, BP, LP, NOTCH, AP};

  float frequency;
  float q;
  FMODE mode;
  int samplerate;
  
  float a1, a2;
  float b0, b1, b2;
  float k;  
  
  float del1, del2; 
  
  canonical_sos_filter ();
  canonical_sos_filter (FMODE mode);
  canonical_sos_filter (float frequency, float q, int samplerate, FMODE mode);

  // update filter parameters
  void update_internals ();
    
  // calculate filtered sample
  float calculate (float sample);  
};


/*
 * simple mean filter, to shave the edge off a little ...
 */
struct simple_mean_filter : public filter {
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

  void update_internals();  
  float calculate (float sample);
};



struct nonlinear_filter : public filter {
  float kp, kn, gp, gn;
  float alpha_mix;
  float gain_sc;
  float g_pre, g_post;
  float lp_freq;
  
  canonical_sos_filter lp;

  nonlinear_filter(float samplerate);
  
  float nonlinear_transfer(float sample);

  float calculate(float sample);
  void update_internals();
};


/*
 * Filterbanks, mostly for "manual mode..."
 */
struct mean_filterbank {
  int channels;
  simple_mean_filter* filters;
  
  mean_filterbank (int channels, int points);

  ~mean_filterbank();
  
  void update(int points);
  
  void process (int channel, float& sample);
};


// a simple filterbank consisting of several sos filters
struct filterbank {
  
  int bands;
  int channels;
  canonical_sos_filter* fbank;
  bool* fmask;

  filterbank(int channels, int samplerate, float lowcut, float hicut, int bands);

  ~filterbank();
  
  void process (int channel, float& sample);

  void toggle_band (int band);

  void print_bands ();
};

#endif
