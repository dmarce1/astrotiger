# Copyright (c) 2015 Thomas Heller
# Copyright (c) 2015 Dominic Marcello
# Copyright (c) 2018-2019 Parsa Amini
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.12.4)
cmake_policy(VERSION 3.12.4)

# Libraries linked via full path no longer produce linker search paths.
cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0074 NEW)

project(astrotiger CXX C)

# Search path for CMake modules to be loaded by include() and find_package()
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(HPX REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
find_package(Silo REQUIRED)
find_package(Boost REQUIRED)


set(CMAKE_CXX_FLAGS "-Wno-pointer-arith -isystem /home/dmarce1/local/boost/include -DBOOST_ALL_NO_LIB -DBOOST_DISABLE_ASSERTS -O3 -std=c++17 -DNDEBUG -ffast-math -march=native")
#set(CMAKE_CXX_FLAGS "-isystem /home/dmarce1/local/boost/include -DBOOST_ALL_NO_LIB -DBOOST_DISABLE_ASSERTS -g -std=c++17 -ffast-math -march=native")

set(astrotiger_source_files
	src/boxes.cpp
	src/cosmos.cpp
	src/geometry.cpp
	src/gravity.cpp
	src/fileio.cpp
	src/hpx.cpp
	src/hydro_grid.cpp
	src/main.cpp
	src/particles.cpp
	src/polytrope.cpp
	src/output.cpp
	src/options.cpp
	src/rand.cpp
	src/tree.cpp
	src/tree_client.cpp
	src/levels.cpp
)


set(astrotiger_header_files
 	astrotiger/boxes.hpp	
 	astrotiger/channel.hpp	
 	astrotiger/cosmos.hpp	
 	astrotiger/defs.hpp	
 	astrotiger/fileio.hpp	
 	astrotiger/geometry.hpp	
 	astrotiger/gravity.hpp	
 	astrotiger/hpx.hpp	
 	astrotiger/hydro_flux.hpp	
 	astrotiger/hydro_grid.hpp	
 	astrotiger/levels.hpp	
 	astrotiger/multi_array.hpp	
 	astrotiger/output.hpp	
 	astrotiger/options.hpp	
 	astrotiger/particles.hpp	
 	astrotiger/polytrope.hpp	
 	astrotiger/rand.hpp
 	astrotiger/range.hpp
 	astrotiger/tree.hpp
 	astrotiger/tree_client.hpp
 	astrotiger/vect.hpp	
)


add_library(hpx_lite 
	hpx_lite/src/hpx.cpp
	hpx_lite/src/id_type.cpp
	hpx_lite/src/thread.cpp
	hpx_lite/src/future.cpp
	hpx_lite/src/mutex.cpp
	hpx_lite/src/serialize.cpp
)

add_executable( astrotiger_lite
    ${astrotiger_source_files}
)

add_hpx_executable(
  astrotiger
  DEPENDENCIES
    Silo::silo Boost::boost
  SOURCES
    ${astrotiger_source_files}
  HEADERS
    ${astrotiger_header_files}
 )



target_link_libraries (astrotiger_lite /home/dmarce1/local/boost/lib/libboost_program_options.so tcmalloc hpx_lite siloh5 hwloc dl z hdf5)


target_include_directories(hpx_lite PUBLIC ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/hpx_lite)
target_include_directories(astrotiger PUBLIC ${PROJECT_SOURCE_DIR})
target_compile_options(astrotiger PUBLIC "-DUSE_HPX")
target_include_directories(astrotiger_lite PUBLIC ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/hpx_lite)
target_compile_options(astrotiger_lite PUBLIC "-DHPX_LITE")


