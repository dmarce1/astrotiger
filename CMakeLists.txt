cmake_minimum_required(VERSION 3.18)

project(astrotiger C CXX CUDA)

set(ALL_FLAGS -std=c++20 -Wall -Wno-comma-subscript -I${PROJECT_SOURCE_DIR}/include)
set(CXX_FLAGS  ${ALL_FLAGS} -march=native)
set(CUDA_FLAGS ${ALL_FLAGS} -Xcompiler -march=native)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(REQUIRED)

include_directories(/usr/include/hdf5/serial)
link_directories(/usr/lib/x86_64-linux-gnu/hdf5/serial)

set(astrotiger_src
  src/main.cpp
  src/Util.cpp
  src/Radiation.cpp
)

set(astrotiger_header
  include/Array.hpp
  include/Basis.hpp
  include/Constants.hpp
  include/ContainerArithmetic.hpp
  include/Definitions.hpp
  include/EulerState.hpp
  include/Hdf5.hpp
  include/HyperGrid.hpp
  include/Legendre.hpp
  include/Matrix.hpp
  include/MultiIndex.hpp
  include/Octant.hpp
  include/Options.hpp
  include/Quadrature.hpp
  include/Range.hpp
  include/Real.hpp
  include/RungeKutta.hpp
  include/Util.hpp
  include/Zobrist.hpp
)

add_executable(codegen
  src/codegen.cpp
  src/permutation.cpp
  src/Util.cpp
  include/Definitions.hpp
  include/Indent.hpp
  include/Util.hpp
)

target_include_directories(codegen PRIVATE
  ${PROJECT_SOURCE_DIR}/include
)

target_compile_options(codegen PUBLIC
   $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -O0 -g -fno-omit-frame-pointer>>
   $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -O2>>
)

target_link_options(codegen PUBLIC
   $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -rdynamic>>
   $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native>>
   $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:-pg ${CXX_FLAGS} -march=native>>
)

set(GENERATED_SOURCE_DIR
    ${CMAKE_BINARY_DIR}/generated_source
)

file(MAKE_DIRECTORY ${GENERATED_SOURCE_DIR})

set(GENERATED_STAMP "${GENERATED_SOURCE_DIR}/.codegen.stamp")

add_custom_command(
  OUTPUT ${GENERATED_STAMP}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${GENERATED_SOURCE_DIR}
  COMMAND codegen --output-dir ${GENERATED_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} -E touch ${GENERATED_STAMP}
  DEPENDS codegen
  COMMENT "Running codegen → ${GENERATED_SOURCE_DIR}"
  VERBATIM
)

add_custom_target(generate_code
  DEPENDS ${GENERATED_STAMP}
)

file(GLOB_RECURSE GENERATED_SOURCES
  CONFIGURE_DEPENDS
  "${GENERATED_SOURCE_DIR}/*.cpp"
)

add_executable(
   astrotiger
      ${astrotiger_src}
      ${GENERATED_SOURCES}
      ${astrotiger_header}
)


target_link_libraries(astrotiger PUBLIC hdf5)

add_dependencies(astrotiger generate_code)

target_include_directories(astrotiger
   PUBLIC
      ${PROJECT_SOURCE_DIR}/include
      ${GENERATED_SOURCE_DIR}
)

set_property(TARGET astrotiger PROPERTY CUDA_ARCHITECTURES native)
set_property(TARGET astrotiger PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_compile_options(astrotiger PUBLIC
   $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O3>>
   $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O0 -g -fno-omit-frame-pointer>>
   $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -g -pg -DNDEBUG>>
   $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -G>>
   $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} --use_fast_math -DNDEBUG Xcompiler=-Ofast>>
   $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -Xcompiler=-g -pg -DNDEBUG>>
)

target_link_options(astrotiger PUBLIC
   $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -rdynamic>>
)


