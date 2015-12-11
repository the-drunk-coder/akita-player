#include <iostream>
#include <sndfile.hh>
#include "RtAudio.h"

#define CHUNKSIZE 1024

enum RWTYPES {
  INT8,
  INT16,
  INT24,
  INT32,
  FLOAT32,
  FLOAT64
};

template<typename READ_TYPE>
struct file_container {
  READ_TYPE* file_buffer;
  READ_TYPE* chunk_buffer;
  long int offset;
  short channels;
};

template<typename READ_TYPE>
struct params_to_abuse {

  // the file to use !
  file_container<READ_TYPE>* sample_file;

  // the params !
  float frame_increment;
  float sample_repeat;
};

// the parameterized callback function ...
template<typename READ_TYPE, typename WRITE_TYPE>
int abusive_play( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
         double streamTime, RtAudioStreamStatus status, void *userData )
{
  // get the parameter container from the user data ...
  params_to_abuse<READ_TYPE>* params = reinterpret_cast<params_to_abuse<READ_TYPE>*>(userData);
	file_container<READ_TYPE>* file_container = params->sample_file;

  if ( status ) { std::cout << "Stream underflow detected!" << std::endl; }

  // transfer samples from file buffer to output !
  for(int i = 0; i < nBufferFrames * 2; i++){
    ((WRITE_TYPE*) outputBuffer)[i] = file_container->file_buffer[file_container->offset + i];
  }

  // increament read offest
  file_container->offset += nBufferFrames * 2;

  return 0;
}

int main () {
	const char * fname = "test.wav" ;

	SndfileHandle file = SndfileHandle (fname) ;

  std::cout << "Reading file: " << fname << std::endl;
  std::cout << "File format: " << file.format() << std::endl;
  std::cout << "PCM 16 BIT: " << (SF_FORMAT_WAV | SF_FORMAT_PCM_16) << std::endl;
  std::cout << "Frames in file: " << file.frames() << std::endl;
  std::cout << "Samplerate " << file.samplerate() << std::endl;
  std::cout << "Channels: " << file.channels() << std::endl;

  params_to_abuse<int16_t>* params = new params_to_abuse<int16_t>;
  params->sample_file = new file_container<int16_t>;

  params->sample_file->file_buffer = new int16_t[file.frames() * file.channels()];

  int frame_chunksize = CHUNKSIZE * file.channels();
  params->sample_file->chunk_buffer = new int16_t[frame_chunksize];
  params->sample_file->offset = 0;
  params->sample_file->channels = file.channels();

  std::cout << "Chunk Buffer allocated !" << std::endl;
  // Read file into buffer
  long int read_offset = 0;
  while(file.readf(params->sample_file->chunk_buffer, CHUNKSIZE) == CHUNKSIZE ){
    for (int i = 0; i < frame_chunksize; i++){
      params->sample_file->file_buffer[read_offset + i] = params->sample_file->chunk_buffer[i];
    }
    read_offset += frame_chunksize;
  }

  std::cout << "File read!" << std::endl;

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
    dac.openStream( &parameters, NULL, RTAUDIO_SINT16,
                    sampleRate, &bufferFrames, &abusive_play<int16_t, int16_t>, (void*) params);
    dac.startStream();
  }
  catch ( RtAudioError& e ) {
    e.printMessage();
    return 0;
  }

  char input;
  std::cout << "\nPlaying ... press <enter> to quit.\n";
  std::cin.get( input );
  try {
    // Stop the stream
    dac.stopStream();
  }
  catch (RtAudioError& e) {
    e.printMessage();
  }
  if ( dac.isStreamOpen() ) dac.closeStream();

  return 0 ;

}
