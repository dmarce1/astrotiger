cmake_minimum_required(VERSION 3.18)

project(astrotiger LANGUAGES CXX CUDA)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-std=c++20)
endif()

# ------------------------------------------------------------------------------
# FLAGS
# ------------------------------------------------------------------------------
set(ALL_FLAGS -std=c++23 -Wall -Wno-comma-subscript -I${PROJECT_SOURCE_DIR}/include)
set(CXX_FLAGS  ${ALL_FLAGS} -fext-numeric-literals -march=native)
set(CUDA_FLAGS ${ALL_FLAGS} -Xcompiler -march=native)

# ------------------------------------------------------------------------------
# DEPENDENCIES
# ------------------------------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(HPX  REQUIRED)


# ------------------------------------------------------------------------------
# SOURCES
# ------------------------------------------------------------------------------
set(astrotiger_src
  src/main.cpp
  src/Constants.cpp
  src/FermiDirac.cpp
  src/Options.cpp
  src/Util.cpp
  src/Radiation.cpp
  src/SodShockTube.cpp
  src/Symbolic.cpp
)

set(astrotiger_header
   include/AutoDiff.hpp		 
   include/Constants.hpp		  
   include/DifferentialPolynomial.hpp  
   include/Face.hpp	       
   include/FwdAutoDiff.hpp 
   include/Hdf5.hpp		 
   include/Indent.hpp	 
   include/Multidices.hpp    
   include/Octogrid.hpp	   
   include/Polynomial.hpp      
   include/Range.hpp     
   include/RungeKutta.hpp   
   include/TriIndex.hpp
   include/BellPolynomial.hpp	
   include/ContainerArithmetic.hpp  
   include/DoubleReal.hpp	     
   include/FermiDirac.hpp   
   include/Hypergeometric.hpp 
   include/Integrate.hpp  
   include/MultiIndex.hpp    
   include/Options.hpp	 
   include/Quadrature.hpp    
   include/Rational.hpp 
   include/stateCodegen.hpp 
   include/Units.hpp
   include/BiconjugateGradient.hpp  
   include/Definitions.hpp	  
   include/EulerState.hpp	    
   include/HalfInteger.hpp 
   include/HyperSubgrid.hpp  
   include/Matrix.hpp	 
   include/MultiPrecision.hpp
   include/PadeApproximant.hpp 
   include/RadiationState.hpp 
   include/Radiation.hpp 
   include/Real.hpp    
   include/Symbolic.hpp   
   include/Util.hpp
)


# ------------------------------------------------------------------------------
# CODEGEN
# ------------------------------------------------------------------------------
add_executable(codegen
  	src/codegen.cpp
  	src/stateCodegen.cpp
  	src/Util.cpp
  	include/Definitions.hpp
  	include/Indent.hpp
  	include/stateCodegen.hpp
	include/Symbolic.hpp
 	include/Util.hpp
)

target_include_directories(codegen PRIVATE
  ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(codegen PUBLIC
    gmp
	stdc++_libbacktrace
    quadmath    
    symengine
)

target_compile_options(codegen PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -O0 -g -fno-omit-frame-pointer>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -O2>>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -g -pg -march=native -O2 -DNDEBUG>>
)

target_link_options(codegen PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -rdynamic>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native>>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -pg -march=native>>
)

# ------------------------------------------------------------------------------
# GENERATED CODE
# ------------------------------------------------------------------------------
set(GENERATED_SOURCE_DIR ${CMAKE_BINARY_DIR}/generated_source)
file(MAKE_DIRECTORY ${GENERATED_SOURCE_DIR})
set(GENERATED_DG_TRANSFORMS "${GENERATED_SOURCE_DIR}/dgTransforms.hpp")

add_custom_command(
  OUTPUT  "${GENERATED_DG_TRANSFORMS}"
          "${GENERATED_SOURCE_DIR}/.codegen.stamp"
  COMMAND "${CMAKE_COMMAND}" -E make_directory "${GENERATED_SOURCE_DIR}"
  COMMAND codegen --output-dir "${GENERATED_SOURCE_DIR}"
  COMMAND "${CMAKE_COMMAND}" -E compare_files
          "${GENERATED_DG_TRANSFORMS}" "${GENERATED_DG_TRANSFORMS}"
  COMMAND "${CMAKE_COMMAND}" -E touch
          "${GENERATED_SOURCE_DIR}/.codegen.stamp"
  DEPENDS codegen
  COMMENT "Running codegen → ${GENERATED_SOURCE_DIR}"
  VERBATIM
)

add_custom_target(generate_code
  DEPENDS "${GENERATED_DG_TRANSFORMS}"
)

file(GLOB_RECURSE GENERATED_SOURCES
  CONFIGURE_DEPENDS
  "${GENERATED_SOURCE_DIR}/*.cpp"
)

# ------------------------------------------------------------------------------
# ASTROTIGER EXECUTABLE
# ------------------------------------------------------------------------------
add_hpx_executable(
  astrotiger
  DEPENDENCIES
    HDF5::HDF5
    stdc++_libbacktrace
    quadmath
    mpfr
    gmp  
  SOURCES
    ${astrotiger_src}
    ${GENERATED_SOURCES}
  HEADERS
    ${astrotiger_header}
)

add_dependencies(astrotiger generate_code)

target_include_directories(astrotiger
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${GENERATED_SOURCE_DIR}
)

set_property(TARGET astrotiger PROPERTY CUDA_ARCHITECTURES native)
set_property(TARGET astrotiger PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_compile_options(astrotiger PUBLIC
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O3 -g -fopt-info-vec-optimized >>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O0 -g -fno-omit-frame-pointer>>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -g -pg -O2 -DNDEBUG>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -G>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} --use_fast_math -DNDEBUG Xcompiler=-Ofast>>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -Xcompiler=-g -pg -DNDEBUG>>
)

target_link_options(astrotiger PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -rdynamic>>
)

# ------------------------------------------------------------------------------
# TESTS
# ------------------------------------------------------------------------------
# Add Catch2
add_subdirectory(extern/Catch2)

# Define test sources
set(test_sources
  include/FwdAutoDiff.hpp 
  include/Matrix.hpp 
  include/Radiation.hpp 
  include/Units.hpp 
  tests/test_con2prim.cpp
  tests/test_implicit_radiation.cpp
  src/Radiation.cpp
  src/Util.cpp
)

add_executable(unit_tests ${test_sources} )

target_link_libraries(unit_tests
  PRIVATE
    Catch2::Catch2WithMain
    HDF5::HDF5
    stdc++_libbacktrace
    quadmath
    mpfr
    gmp
)

target_include_directories(unit_tests
  PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${GENERATED_SOURCE_DIR}
)

target_compile_options(unit_tests PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O0 -g>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O3>>
)

include(CTest)
include(Catch)
catch_discover_tests(unit_tests)

