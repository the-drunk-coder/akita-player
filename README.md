# akita-player
Stupid little program that you can use to create noise abusing low-level audio parameters ...

Needs the following libraries:
* [**libsndfile**](https://github.com/erikd/libsndfile)
* [**boost::program_options**](http://www.boost.org/doc/libs/1_61_0/doc/html/program_options/tutorial.html)
* [**boost::lockfree**](http://www.boost.org/doc/libs/1_61_0/doc/html/lockfree.html)
* [**jack**](http://www.jackaudio.org)

It'll currently only work with JACK, unless you want to use RAW mode all the time.

If you want to compile it with another API, please refer to the RtAudio documentation!

Assuming you have all libraries available, you can build *akita* with **cmake**:

```
akita/$ mkdir build && cd build
akita/build/$ cmake ..
akita/build/$ make
akita/build/$ make install
```

Can be used like an audio player (for WAV or FLAC), but with a little twist, as you can (ab-)use incompatible read/write/stream-types, underfilled buffers, down- and resampling and the like ... just toy around with it a little (but at **low volumes**, as the outcome might be **ear-shattering**...).

Some options or combinations might work all the time ... 

For notes on usage, type s.th. like once you've compiled it !

```./akita --help```

This software is public domain !
