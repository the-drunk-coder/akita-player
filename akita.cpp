#include <iostream>
#include <algorithm>
#include <iterator>
#include <boost/program_options.hpp>
#include <mutex>              
#include <condition_variable>
#include <string>
#include <functional>
#include <climits>
#include <sys/time.h>
#include <iomanip>
#include <lo/lo.h>
#include <lo/lo_cpp.h>
#include "getch.h"
#include "RtAudio.h"
#include "akita_structures.h"
#include "akita_actions.h"

namespace po = boost::program_options;

// mutex and cv to block audio thread 
std::mutex mtx;
std::condition_variable cv;

/*
 * RtAudio callback functions for source- and filter clients.
 */

// The parameterized generator callback function ...
template <typename READ_TYPE, typename WRITE_TYPE>
int source_callback(void *outputBuffer, void *inputBuffer,
                 unsigned int nBufferFrames, double streamTime,
                 RtAudioStreamStatus status, void *userData) {
  
  // get the parameter container from the user data ...
  source_params<READ_TYPE> *spar = reinterpret_cast<source_params<READ_TYPE>*>(userData);  
  file_container<READ_TYPE>& fc = spar->fc;

  WRITE_TYPE *out_buf = (WRITE_TYPE *) outputBuffer;

  if (status) { std::cout << "Stream underflow detected!" << std::endl; }
  
  uint32_t block_pos = 0;
  if (spar->iface == OSC) {    
    if (spar->current_event->state != akita_play_event::IN_PROGRESS) {      
      for (int i = 0; i < nBufferFrames * fc.channels; i++) {		
	if (spar->current_event->state == akita_play_event::NEW) {
	    fc.update_range(*spar->current_event);
	    spar->offset = fc.start_sample;
	    spar->state = PLAY_EVENT;
	    spar->current_event->state = akita_play_event::IN_PROGRESS;	    
	    break;
	}
	out_buf[block_pos++] = 0;			
      }      
    }
  }
                       
  // block source thread ... might have an interesting effect ...
  if(spar->state == BLOCK || spar->state == LOOP_BLOCK){
    std::unique_lock<std::mutex> lck(mtx);
    cv.wait(lck);
    return 0;
  }
  
  // SILENCE ! I KILL YOU !!
  if (spar->state == SILENCE || spar->state == LOOP_SILENCE) {
    for (int i = 0; i < nBufferFrames * spar->buffer_cut; i++) {
      *out_buf++ = 0;
    }
    // in this case, offset won't be modified
    return 0;
  }

  // transfer samples from file buffer to output !
  float file_buf_pos = 0;
  for (uint32_t i = block_pos; i < nBufferFrames * fc.channels; i++) {   
    // loop in case we hit the file's end
    if (spar->offset + (int) file_buf_pos > fc.end_sample) {      
      spar->offset = fc.start_sample;
      file_buf_pos = 0;
    }

    if (spar->state == PLAY_EVENT && spar->current_event->samples_played >= spar->current_event->length_samples ) {
      spar->state = SILENCE;
      spar->current_event->state = akita_play_event::FINISHED;
      break;
    }
    
    // copy samples        
    
    out_buf[i] = fc.file_buffer[spar->offset + (int) file_buf_pos];    
    file_buf_pos += 1.0 / spar->samplerate_mod;
    
    if (spar->state == PLAY_EVENT){
      spar->current_event->samples_played++;
    } 
  }

  // sample repetition
  if (spar->sample_repeat > 1) {
    for (uint32_t i = block_pos; i < (nBufferFrames * fc.channels) / spar->sample_repeat; i++) {  
      for(int j = 0; j < spar->sample_repeat; j++){
	out_buf[i+j] = out_buf[i];
      }
    }
  } 
  
  // buffercutting
  for (uint32_t i = nBufferFrames * spar->buffer_cut; i < nBufferFrames * fc.channels; i++) {           
    out_buf[i] = 0.0;    
  }
  
  // random sample shootout
  if(spar->fuzziness > 0.0){
    for (uint32_t i = block_pos; i < nBufferFrames * fc.channels; i++) {   
      if(rand() / (float) INT_MAX < spar->fuzziness){
	out_buf[i] = 0.0;    
      }
    }
  }
  
  // increament read offset
  spar->offset += ((float) (nBufferFrames - block_pos) * spar->offset_cut) / spar->sample_repeat;
  if((spar->state == LOOP) && (spar->offset >= spar->loop_end)){
    spar->offset = spar->loop_start.load();
  }
   
  return 0;
}

// the callback function for the filter thread
int filter_callback(void *outputBuffer, void *inputBuffer,
                 unsigned int nBufferFrames, double streamTime,
	   RtAudioStreamStatus status, void *userData) {

  filter_params *fpar = reinterpret_cast<filter_params*>(userData);
    
  float* in_buf = (float *) inputBuffer;
  float* out_buf = (float *) outputBuffer;

  long int frame_offset = 0;
  for (int i = 0; i < nBufferFrames; i++) {
    float current_sample = in_buf[i];

    if(fpar->flippiness > 0.0) {
      current_sample = random_flip(current_sample, fpar->flippiness);
    }
   
    if(fpar->mode == PMODE::MILD){
      // weed out impossible values ... while this kinda takes some of the edge off,
      // it enables us to use some stuff like filters and gain ...
      if (isnan(current_sample)) current_sample = 0;    
      if (current_sample < -1.0) current_sample = -1.0;
      if (current_sample > 1.0) current_sample = 1.0;
    }  
    
    if(fpar->mean_filter_on){
      fpar->m_fbank.apply(0, current_sample);
    }
    
    if (fpar->filterbank_on) {    
      fpar->fbank.apply(0, current_sample);      
    }

    if (fpar->lowpass_on) {    
      fpar->lowpass.process(current_sample);      
    }
    
    if(fpar->mode == PMODE::MILD){
      current_sample *= fpar->gain;
    }
    
    fpar->frame_buffer[fpar->pan_offset] = current_sample * (1.0 - fpar->pan_ratio);
    fpar->frame_buffer[fpar->pan_offset + 1] = current_sample * fpar->pan_ratio;

    if (fpar->reverb_on) {      
      fpar->rev.tick(fpar->frame_buffer[fpar->pan_offset], fpar->frame_buffer[fpar->pan_offset + 1]);
      fpar->frame_buffer[fpar->pan_offset] = fpar->rev.lastOut(0);
      fpar->frame_buffer[fpar->pan_offset + 1] = fpar->rev.lastOut(1);
    }

    for (int j = 0; j < fpar->channels; j++) {
      out_buf[frame_offset + j] = fpar->frame_buffer[j];
    }
    frame_offset += fpar->channels;    
  }
  
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

// initialize the command line parameter parser
po::options_description init_opts(int ac, char *av[], po::variables_map& vm,
                                  options_container& opts) {

  po::options_description desc("Parameters to use in a creative way");
  desc.add_options()
    ("help", "Display this help!")
    ("interface", po::value<INTERFACE>(&opts.iface)->default_value(PLAIN), "Interaction mode")
    ("input-file", po::value<std::string>(&opts.filename)->default_value(""), "The input file - WAV of FLAC!")    
    ("init-state", po::value<PSTATE>(&opts.initial_state)->default_value(PLAY), "Initial state!")    
    ("init-mode", po::value<PMODE>(&opts.initial_mode)->default_value(MILD), "Initial mode!")
    ("init-gain", po::value<float>(&opts.initial_gain)->default_value(0.5), "Initial gain (default: 0.5)!")
    ("start", po::value<float>(&opts.start)->default_value(0.0), "Starting point within the sample, relative to length!")
    ("end", po::value<float>(&opts.end)->default_value(1.0), "End point within the sample, relative to length!")
    ("samplerate-mod", po::value<float>(&opts.samplerate_mod)->default_value(1), "Modifiy samplerate!")
    ("sample-repeat", po::value<int>(&opts.sample_repeat)->default_value(1), "Repeat every sample n times!")
    ("buffer-cut", po::value<float>(&opts.buffer_cut)->default_value(2), "Don't fill source output buffer completely!")
    ("offset-mod", po::value<float>(&opts.offset_cut)->default_value(2), "Modify offset increment (chunk size read from buffer)!")
    ("fuzziness", po::value<float>(&opts.fuzziness)->default_value(0.0), "Create fuzziness by removing random samples with a certain probability!")
    ("flippiness", po::value<float>(&opts.flip_prob)->default_value(0.0), "Create different fuzziness by flipping bits a certain probability!")
    ("read-type", po::value<RWTYPES>(&opts.read_type)->default_value(SHORT), "Type used to read from audio file!")
    ("write-type", po::value<RWTYPES>(&opts.write_type)->default_value(SHORT), "Type used to write to audio buffer!")
    ("stream-type", po::value<RWTYPES>(&opts.stream_type)->default_value(SHORT), "Type used for the audio stream!")
    ("pan", po::value<float>(&opts.pan)->default_value(0.5), "pan!")
    ("out-channels", po::value<int>(&opts.out_channels)->default_value(2), "Number of output channels!")
    ("mean-filter", po::value<int>(&opts.mean_filter_points)->default_value(0), "Apply mean filter to shave the edge off a little!")
    //("mono", po::value<bool>(&opts.mono)->default_value(true), "Mixdown to mono!")
    ("reverb", po::value<float>(&opts.reverb_mix)->default_value(0.4), "Reverb level!")
    ("udp-port", po::value<int>(&opts.udp_port)->default_value(19456), "UDP Port for OSC mode!")    
    ;
  // ----- end options ... what kind of syntax is this ??

  po::positional_options_description p;
  p.add("input-file", -1);

  po::store(po::command_line_parser(ac, av).options(desc).positional(p).run(),vm);
  po::notify(vm);

  return desc;
}

/*
 * Function to stop audio
 */ 
void stop_audio(RtAudio& source, RtAudio& filter, bool raw) {
  source.stopStream();
  source.closeStream();

  if (!raw){
    filter.stopStream();
    filter.closeStream();
  }
}

template<typename READ_TYPE>
void handle_osc_input(source_params<READ_TYPE>& spar, filter_params& fpar, options_container& opts){
  
  // init osc server
  lo::ServerThread st(opts.udp_port);

  // an exception might be nice here ...
  if (!st.is_valid()) {
    std::cout << "Can't start udp server." << std::endl;
    return;
  }

  // message handlers 
  st.add_method("/akita/play", "fiffiffffiff",
		[&spar, &fpar, &opts](lo_arg **argv, int count) {
		  akita_actions::change_gain(spar, fpar, argv[2]->f);
		  akita_actions::change_reverb_mix(spar, fpar, argv[3]->f);
		  fpar.mean_filter_on = argv[4]->i;
		  akita_actions::change_lowpass(spar, fpar, argv[5]->f, argv[6]->f);
		  akita_actions::change_flippiness(spar, fpar, argv[7]->f);
		  akita_actions::change_fuzziness(spar, fpar, argv[8]->f);		  
		  akita_actions::change_pan(spar, fpar, argv[10]->f);

		  // otherwise, the event is still in progress
		  if (!(spar.current_event != NULL && spar.current_event->state != akita_play_event::FINISHED)){
		    std::cout << "A K I T A - instance at " << opts.udp_port << " RECIEVED play EVENT ! " << std::endl;		    
		    spar.current_event.reset(new akita_play_event(argv[0]->f, argv[1]->i));
		    
		    akita_actions::change_samplerate_mod(spar, fpar, argv[11]->f);
		    akita_actions::change_sample_repeat(spar, fpar, argv[9]->i);
		    return 0;
		  } else {
		    std::cerr << "A K I T A - IGNORED play EVENT !" << std::endl;
		    return 0;
		  } 		  		  
		});

  // message handlers 
  st.add_method("/akita/param", "ffiffffiff",
		[&spar, &fpar](lo_arg **argv, int count) {

		  std::cout << "A K I T A - RECIEVED param EVENT !" << std::endl;
		  
		  akita_actions::change_gain(spar, fpar, argv[0]->f);
		  akita_actions::change_reverb_mix(spar, fpar, argv[1]->f);
		  fpar.mean_filter_on = argv[2]->i;
		  akita_actions::change_lowpass(spar, fpar, argv[3]->f, argv[4]->f);
		  akita_actions::change_flippiness(spar, fpar, argv[5]->f);
		  akita_actions::change_fuzziness(spar, fpar, argv[6]->f);
		  akita_actions::change_sample_repeat(spar, fpar, argv[7]->i);
		  akita_actions::change_pan(spar, fpar, argv[8]->f);
		  akita_actions::change_samplerate_mod(spar, fpar, argv[9]->f); 		  		   		  		  
		});
  
  // message handlers 
  st.add_method("/akita/play", "fi",
		[&spar](lo_arg **argv, int count) {
		  // otherwise, the event is still in progress
		  if (!(spar.current_event != NULL && spar.current_event->state != akita_play_event::FINISHED)){
		    std::cout << "RECIEVED EVENT !" << std::endl;		    
		    spar.current_event.reset(new akita_play_event(argv[0]->f, argv[1]->i));
		    return 0;
		  } else {
		    std::cout << "IGNORED EVENT !" << std::endl;
		    return 0;
		  } 		  		  
		});
  
  // message handlers 
  st.add_method("/akita/param/reverb", "f",
		[&spar, &fpar](lo_arg **argv, int count) {		  		  
		  akita_actions::change_reverb_mix(spar, fpar, argv[0]->f);
		});

  // message handlers 
  st.add_method("/akita/param/mean_filter", "i",
		[&spar, &fpar](lo_arg **argv, int count) {
		  fpar.mean_filter_on = argv[0]->i;
		});

  // message handlers 
  st.add_method("/akita/param/lowpass", "ff",
		[&spar, &fpar](lo_arg **argv, int count) {		  
		  akita_actions::change_lowpass(spar, fpar, argv[0]->f, argv[1]->f);
		});

  // message handlers 
  st.add_method("/akita/param/flippiness", "f",
		[&spar, &fpar](lo_arg **argv, int count) {
		  akita_actions::change_flippiness(spar, fpar, argv[0]->f);
		});

  // message handlers 
  st.add_method("/akita/param/fuzziness", "f",
		[&spar, &fpar](lo_arg **argv, int count) {
		  akita_actions::change_fuzziness(spar, fpar, argv[0]->f);
		});


  // message handlers 
  st.add_method("/akita/param/gain", "f",
		[&spar, &fpar](lo_arg **argv, int count) {
		  akita_actions::change_gain(spar, fpar, argv[0]->f);
		});


  // mutex and cv to block audio thread 
  std::mutex osc_mtx;
  std::condition_variable osc_cv;
  
  // message handlers 
  st.add_method("/akita/quit", "",
		[&osc_cv](lo_arg **argv, int count) {
		  std::cout << "Quitting akita instance!" << std::endl;
		  osc_cv.notify_all();
		});
  
  // start server 
  st.start();

  // wait for quit signal
  std::unique_lock<std::mutex> lck(osc_mtx);
  osc_cv.wait(lck);    
}

template<typename READ_TYPE>
void handle_keyboard_input(source_params<READ_TYPE>& spar, filter_params& fpar ){
  char input;
  // main loop
  
  while((input = getch()) != 'q'){       
    switch(input) {    
      // filter control
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
    case '0':
      akita_actions::toggle_filter_band(spar, fpar, input);      
      break;
    case 'f':      
      akita_actions::toggle_filterbank(spar, fpar);
      break;
    case 'r':
      akita_actions::toggle_reverb(spar, fpar);
      break;
    case 'g':
      akita_actions::toggle_mean_filter(spar, fpar);
      break;
    // gain control
    case 'd':           
      akita_actions::change_gain(spar, fpar, fpar.gain + 0.05);
      break;
    case 'c':      
      akita_actions::change_gain(spar, fpar, fpar.gain - 0.05);
      break;
      // fuzziness control
    case 'a':           
      akita_actions::change_fuzziness(spar, fpar, spar.fuzziness + 0.01);
      break;      
    case 'y':      
      akita_actions::change_fuzziness(spar, fpar, spar.fuzziness - 0.01);
      break;
      // samplerate control
    case 's':           
      akita_actions::change_sample_repeat(spar, fpar, spar.sample_repeat + 1);
      break;
    case 'x':      
      akita_actions::change_sample_repeat(spar, fpar, spar.sample_repeat - 1);
      break;
      // thread blocking control
    case 'm':      
      akita_actions::toggle_mute(spar, fpar);
      break;    
    case ' ':      
      akita_actions::toggle_loop_state(spar, fpar);
      break;
    case 'b':      
      akita_actions::toggle_block(spar, fpar, cv);
      break;
    default:
      std::cout << "COMMAND NOT ACCEPTED!" << std::endl;
    }
    
  }
}


/*
* The Audio initializing function !
*/
template <typename READ_TYPE, typename WRITE_TYPE>
int handle_audio(options_container& opts) {

  if (opts.filename == "") {
    std::cout << "\nPlease specify input file!\n";
    return EXIT_FAILURE;
  }
  
  source_params<READ_TYPE> spar(opts);

  if (spar.fc.state != file_container<READ_TYPE>::READY) {
    std::cout << "File \"" << opts.filename << "\" is not valid!" << std::endl;
    return EXIT_FAILURE;
  }

  if(opts.iface == OSC) {
    std::cout << "Listening on OSC port: " << opts.udp_port << std::endl << std::endl;
  }
  
  spar.fc.print_file_info();
  filter_params fpar(opts);
  fpar.rev.setSampleRate(spar.fc.samplerate);
  fpar.rev.setEffectMix(opts.reverb_mix);
  
  
  std::cout << "Current Parameters:" << std::endl;
  std::cout << "  Read type:     " << rwtypes_strings[opts.read_type] << std::endl;
  std::cout << "  Write type:    " << rwtypes_strings[opts.write_type] << std::endl;
  std::cout << "  Stream type:   " << rwtypes_strings[opts.stream_type] << std::endl;
  std::cout << "  Sample repeat: " << opts.sample_repeat << std::endl;
  std::cout << "  Buffer mod:    " << opts.buffer_cut << std::endl;
  std::cout << "  Offset mod:    " << opts.offset_cut << std::endl << std::endl;
  
  // get current pid
  std::ostringstream str_pid;
  str_pid << ::getpid();

  // initialize audio streams
  RtAudio source;
  RtAudio filter;

  if (filter.getDeviceCount() < 1) {
    std::cout << "\nNo audio devices found!\n";
    return EXIT_FAILURE;
  }
  
  RtAudio::StreamParameters rt_source_parameters;
  
  // in raw mode, all the beatiful glitches will go directly to the DAC ... good luck !
  if(opts.initial_mode != PMODE::RAW){
    RtAudio::StreamParameters rt_filter_out_parameters;
    rt_filter_out_parameters.deviceId = filter.getDefaultOutputDevice();
    rt_filter_out_parameters.nChannels = opts.out_channels;
    rt_filter_out_parameters.firstChannel = 0;

    RtAudio::StreamParameters rt_filter_in_parameters;
    rt_filter_in_parameters.deviceId = filter.getDefaultOutputDevice();
    rt_filter_in_parameters.nChannels = spar.fc.channels;
    rt_filter_in_parameters.firstChannel = 0;

    RtAudio::StreamOptions rt_filter_opts;
    // rt_filter_opts.flags = RTAUDIO_NONINTERLEAVED;
    rt_filter_opts.streamName = "akita-filter-";
    rt_filter_opts.streamName.append(str_pid.str());
    rt_filter_opts.autoConnectInput = false;

    // initialize filter client
    try {
      filter.openStream(&rt_filter_out_parameters, &rt_filter_in_parameters, RTAUDIO_FLOAT32, opts.samplerate,
			&opts.buffer_frames, &filter_callback,
			(void *)&fpar, &rt_filter_opts);
      filter.startStream();
      //set timestamp
      //timeval tv;
      //gettimeofday(&tv, 0);
      //double stream_time =  2208988800.0 + (double) tv.tv_sec + (double)tv.tv_usec / 1000000;// + tv.tv_sec;
      //std::cout << "akita filter stream startup time: " << doubleToText(stream_time) << std::endl;
      //filter.setStreamTime(stream_time);
    } catch (RtAudioError &e) {
      e.printMessage();
      return EXIT_FAILURE;
    }
    
    //find filter thread
    for(unsigned int i = 0; i < source.getDeviceCount(); i++){
      if (source.getDeviceInfo(i).name == rt_filter_opts.streamName){
	rt_source_parameters.deviceId = i;
      }
    }
  } else {
    rt_source_parameters.deviceId = source.getDefaultOutputDevice();;
  }

  rt_source_parameters.nChannels = spar.fc.channels;
  rt_source_parameters.firstChannel = 0;

  RtAudio::StreamOptions source_opts;
  source_opts.streamName = "akita-noise-source-";
  source_opts.streamName.append(str_pid.str());

  // initialize noise source
  try {
    source.openStream(&rt_source_parameters, NULL, rwtypes_rtaudio[opts.stream_type], opts.samplerate,
		       &opts.buffer_frames, &source_callback<READ_TYPE, WRITE_TYPE>,
		       (void *)&spar, &source_opts);
    
    source.startStream();
    // set timestamp
    //timeval tv;
    //gettimeofday(&tv, 0);
    //double stream_time =  2208988800.0 + (double) tv.tv_sec + (double)tv.tv_usec / 1000000;// + tv.tv_sec;
    //std::cout << "akita source stream startup time: " << doubleToText(stream_time) << std::endl;
    //source.setStreamTime(stream_time);
  } catch (RtAudioError &e) {
    e.printMessage();
    return EXIT_FAILURE;
  }
    
  if(opts.iface == OSC){
    std::cout << "Listening for OSC input!" << std::endl;
    handle_osc_input(spar, fpar, opts);
  } else {
    std::cout << "Playing ... press q to quit" << std::endl;
    handle_keyboard_input(spar, fpar);
  }

  try {
    stop_audio(source, filter, (opts.initial_mode == PMODE::RAW));
  } catch(RtAudioError &e){
    e.printMessage();
    return EXIT_FAILURE;
  }

  std::cout << "\nQuitting, bye!.\n";
  
  return EXIT_SUCCESS;
}

typedef std::pair<RWTYPES, RWTYPES> init_key;

typedef std::function<int (options_container&)> audio_handler;

std::map<init_key, audio_handler> handlers {
  { init_key(SHORT, SHORT), handle_audio<int16_t, int16_t> },
  { init_key(SHORT, FLOAT), handle_audio<int16_t, float> },
  { init_key(SHORT, DOUBLE), handle_audio<int16_t, double> },
  { init_key(SHORT, UCHAR), handle_audio<int16_t, uint8_t> },
  { init_key(FLOAT, SHORT), handle_audio<float, int16_t> },
  { init_key(FLOAT, FLOAT), handle_audio<float, float> },
  { init_key(FLOAT, DOUBLE), handle_audio<float, double> },
  { init_key(FLOAT, UCHAR), handle_audio<float, uint8_t> },
  { init_key(DOUBLE, SHORT), handle_audio<double, int16_t> },
  { init_key(DOUBLE, FLOAT), handle_audio<double, float> },
  { init_key(DOUBLE, DOUBLE), handle_audio<double, double> },
  { init_key(DOUBLE, UCHAR), handle_audio<double, uint8_t> },
};

/*
 * MAIN ROUTINE!!
 */
int main(int ac, char *av[]) {
  // Salutations!
  std::cout << "\n A K I T A - create noise abusing low-level audio parameters! \n" << std::endl;

  options_container opts;
  
   // initialize command line options
  po::variables_map vm;
  po::options_description desc = init_opts(ac, av, vm, opts);

  // help output
  if (vm.count("help")) {
    std::cout << "usage: akita <file> [options] \n" << std::endl;
    std::cout << desc;
    std::cout << "\nRead/Write/Stream types options:" << std::endl;
    std::cout << "  uchar      8 Bit, unsigned char (not possible as read type!)" << std::endl;
    std::cout << "  short      16 Bit, short (default)" << std::endl;
    std::cout << "  float      32 Bit, float" << std::endl;
    std::cout << "  double     64 Bit, double" << std::endl;
    return EXIT_SUCCESS;
  }

  if (opts.iface == ADVANCED) {
    std::cout << "advanced mode not yet implemented" << std::endl;
    return EXIT_FAILURE;
  } else {
    return handlers[init_key(opts.read_type, opts.write_type)](opts);
  }
}
