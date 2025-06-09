cmake_minimum_required(VERSION 3.18)

project(astrotiger CXX CUDA)

set(
   BUILD_OPTIONS_FILE
      "${CMAKE_SOURCE_DIR}/BuildOptions.cmake"
   CACHE FILEPATH
      "Path to a CMake file that sets preprocessor constants"
)

set(ALL_FLAGS -std=c++23 -Wall -Wno-comma-subscript -I${PROJECT_SOURCE_DIR}/include)
set(CXX_FLAGS  ${ALL_FLAGS} -march=native)
set(CUDA_FLAGS ${ALL_FLAGS} -Xcompiler -march=native)

set(DIMENSION_COUNT 3       CACHE STRING "Number of spatial dimensions (must be >= 1)")
set(SIMD_WIDTH      4       CACHE STRING "Width of the native simd type (must be >= 1)")
set(MODE_COUNT      3       CACHE STRING "Number of modes (must be >= 1)")
set(REAL_DATA_TYPE "double" CACHE STRING "Compile‐time type for all real‐valued data. Default is double.")

if(EXISTS "${BUILD_OPTIONS_FILE}")
    message(STATUS "Including build options from ${BUILD_OPTIONS_FILE}")
    include("${BUILD_OPTIONS_FILE}")
endif()
if(DIMENSION_COUNT LESS 1)
  message(FATAL_ERROR "DIMENSION_COUNT must be at least 1 (got ${DIMENSION_COUNT})")
endif()
if(SIMD_WIDTH LESS 1)
  message(FATAL_ERROR "SIMD_WIDTH must be at least 1 (got ${SIMD_WIDTH})")
endif()
if(MODE_COUNT LESS 1)
  message(FATAL_ERROR "MODE_COUNT must be at least 1 (got ${MODE_COUNT})")
endif()

add_compile_definitions(
   DIMENSION_COUNT=${DIMENSION_COUNT}
   MODE_COUNT=${MODE_COUNT}
   SIMD_WIDTH=${SIMD_WIDTH}
   REAL_DATA_TYPE=${REAL_DATA_TYPE}
)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(HPX  REQUIRED)

set(astrotiger_src
  src/DG.cpp
  src/main.cpp
  src/HyperGrid.cpp
  src/Options.cpp
  src/Quadrature.cpp
  src/TrialSpace.cpp
  src/Util.cpp
)

set(astrotiger_header
  include/Definitions.hpp
  include/HyperGrid.hpp
  include/Options.hpp
  include/Triangular.hpp
  include/TrialSpace.hpp
  include/Quadrature.hpp
  include/Util.hpp
  include/Valarray.hpp
)

add_hpx_executable(
   astrotiger
   DEPENDENCIES
      HDF5::HDF5
      siloh5
      stdc++_libbacktrace
   SOURCES
      ${astrotiger_src}
      ${GENERATED_SOURCES}
   HEADERS
      ${astrotiger_header}
)

set_property(TARGET astrotiger PROPERTY CUDA_ARCHITECTURES native)
set_property(TARGET astrotiger PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_include_directories(astrotiger
   PUBLIC
      ${PROJECT_SOURCE_DIR}/include
)

target_compile_options(astrotiger PUBLIC
   $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O0 -g -fno-omit-frame-pointer>>
   $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O3>>
   $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -g -pg -DNDEBUG >>
   $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -G>>
   $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} --use_fast_math -DNDEBUG Xcompiler=-Ofast>>
   $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -Xcompiler=-g -pg -DNDEBUG>>
)

target_link_options(astrotiger PUBLIC
   $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -rdynamic>>
)

