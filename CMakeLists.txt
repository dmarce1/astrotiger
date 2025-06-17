cmake_minimum_required(VERSION 3.18)
project(astrotiger CXX CUDA)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(HPX REQUIRED)

set(ALL_FLAGS -Wall -Wno-comma-subscript -I${PROJECT_SOURCE_DIR}/include )
set(CXX_FLAGS  ${ALL_FLAGS} -std=c++23 -march=native )
set(CUDA_FLAGS ${ALL_FLAGS} -Xcompiler -march=native )

set(astrotiger_src
	src/main.cpp
	src/Options.cpp
	src/Util.cpp
) 

set(astrotiger_header
	include/Array.hpp
	include/ContainerArithmetic.hpp
	include/EulerState.hpp
	include/Hdf5.hpp
	include/HyperGrid.hpp
	include/HydroGrid.hpp
	include/Legendre.hpp
	include/Matrix.hpp
	include/MultiIndex.hpp
	include/Octant.hpp
	include/Options.hpp
	include/Range.hpp
	include/Real.hpp
	include/RungeKutta.hpp
	include/Util.hpp
	include/Zobrist.hpp
) 

add_executable(codegen
  src/codegen.cpp
  src/Util.cpp
  include/Indent.hpp
  include/Util.hpp
)

target_include_directories(codegen PRIVATE
  ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(codegen PUBLIC
   stdc++_libbacktrace
)

target_compile_options(codegen PUBLIC
   $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -O0 -g -fno-omit-frame-pointer -fsanitize=address>>
   $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -O2>>
)

target_link_options(codegen PUBLIC
   $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -rdynamic -fsanitize=address>>
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

add_hpx_executable(
  astrotiger
  DEPENDENCIES
  	siloh5
	stdc++_libbacktrace
  SOURCES
    ${astrotiger_src}
  HEADERS
    ${astrotiger_header}
)


target_include_directories(astrotiger PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${GENERATED_SOURCE_DIR}
)

set_property(TARGET astrotiger PROPERTY CUDA_ARCHITECTURES native)
set_property(TARGET astrotiger PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_compile_options(
    astrotiger
    PUBLIC
    $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>: ${CXX_FLAGS} -O0 -g -fno-omit-frame-pointer >>
    $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>: ${CXX_FLAGS} -O3 >>
    $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>: ${CXX_FLAGS} -g -pg -DNDEBUG>>
    $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>: ${CUDA_FLAGS} -G >>
    $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CUDA>: ${CUDA_FLAGS} --use_fast_math -DNDEBUG Xcompiler=-Ofast >>
    $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CUDA>: ${CUDA_FLAGS} -Xcompiler=-g -pg -DNDEBUG>>
)

target_link_options(
	astrotiger
  PUBLIC
	$<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>: ${CXX_FLAGS} -rdynamic>>
	$<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>: ${CXX_FLAGS}>>
   	$<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>: -pg ${CXX_FLAGS}>> 
)

target_link_libraries(
	astrotiger 
  PUBLIC 
	HDF5::HDF5
)

target_link_libraries(
	astrotiger 
  PUBLIC
  	$<$<CONFIG:Debug>:stdc++_libbacktrace>
)


add_dependencies(astrotiger generate_code)

