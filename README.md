# akita-player
Stupid little program that you can use to create noise abusing low-level audio parameters ...

Needs [**libsndfile**](https://github.com/erikd/libsndfile), [**boost::program_options**](http://www.boost.org/doc/libs/1_59_0/doc/html/program_options/tutorial.html) and [**RtAudio**](https://github.com/thestk/rtaudio) available!

Assuming you have all libraries available as shared libraries, compile using:

 ```g++ -lboost_program_options -lsndfile -lrtaudio -std=c++11 akita.cpp -o akita ```

Can be used like an audio player (for WAV or FLAC), but with a little twist, as you can (ab-)use incompatible read/write/stream-types, underfilled buffers, down- and resampling and the like ... just toy around with it a little (but at low volumes, as the outcome might be ear-shattering...).

For notes on usage, type s.th. like once you've compiled it !

```./akita --help```

This software is public domain !
