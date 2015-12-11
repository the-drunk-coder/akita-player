#include <iostream>
#include	<sndfile.hh>
#include "RtAudio.h"

#define CHUNKSIZE 1024


struct FileContainer{
  int16_t* file_buffer;
  long int offset;
  short channels;
};

int saw( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
         double streamTime, RtAudioStreamStatus status, void *userData )
{

  int16_t *out_buffer = (int16_t *) outputBuffer;
	FileContainer* file_container = reinterpret_cast<FileContainer*>(userData);

  if ( status ){
    std::cout << "Stream underflow detected!" << std::endl;
  }
  //std::cout << "start!!" << std::endl;
//  std::cout << file_container->file_buffer[file_container->offset + 0] << std::endl;
  // copy from large sample buffer to output buffer!
  for(int i = 0; i < nBufferFrames * 2; i++){

    out_buffer[i] = file_container->file_buffer[file_container->offset + i];
  }
  file_container->offset += nBufferFrames * 2;
//  std::cout << file_container->offset << std::endl;

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

  FileContainer* file_container = new FileContainer;
  file_container->file_buffer = new int16_t[file.frames() * file.channels()];
  std::cout << "Buffer allocated !" << std::endl;
  file_container->offset = 0;
  file_container->channels = file.channels();


  int16_t frame_chunksize = CHUNKSIZE * 2;
  int16_t chunk_buffer[frame_chunksize];
  std::cout << "Chunk Buffer allocated !" << std::endl;
  // Read file into buffer
  long int read_offset = 0;
  while(file.readf(chunk_buffer, CHUNKSIZE) == CHUNKSIZE ){
    for (int i = 0; i < frame_chunksize; i++){
      file_container->file_buffer[read_offset + i] = chunk_buffer[i];
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
                    sampleRate, &bufferFrames, &saw, (void *)file_container);
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
