#cmake_minimum_required (VERSION 3.5)
# use older version for travis compat ...
cmake_minimum_required (VERSION 2.8)

# add project name
project (akita)

# find required libs
find_library(LIBSNDFILE NAMES libsndfile sndfile)
if(NOT LIBSNDFILE)
	message(FATAL_ERROR "Libsndfile not found, please install !")
endif()

# set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -D__UNIX_JACK__")

# internal libraries
add_library(akita_filters STATIC akita_filters.cpp)
add_library(RtAudio STATIC RtAudio.cpp)

# the executable
add_executable(akita akita.cpp)

# link internal libraries
target_link_libraries (akita akita_filters)
target_link_libraries (akita RtAudio)

# link external libraries
target_link_libraries (akita sndfile)
target_link_libraries (akita boost_program_options)
target_link_libraries (akita jack)
target_link_libraries (akita pthread)
target_link_libraries (akita rt)
target_link_libraries (akita stk)

# installation
install(TARGETS akita DESTINATION bin)
