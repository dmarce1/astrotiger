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
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------------------------------------------------------
# DEPENDENCIES
# ------------------------------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(HPX  REQUIRED)

# ------------------------------------------------------------------------------
# IWYU CONFIGURATION
# ------------------------------------------------------------------------------
option(ENABLE_IWYU "Run Include-What-You-Use automatically" OFF)

if(ENABLE_IWYU)
    find_program(IWYU_PATH NAMES include-what-you-use iwyu)
    if(IWYU_PATH)
        message(STATUS "Found IWYU: ${IWYU_PATH}")
        # Default IWYU mapping/behavior
        set(IWYU_ARGS "-Xiwyu;--no_default_mappings")
        set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE "${IWYU_PATH};${IWYU_ARGS}")
    else()
        message(WARNING "IWYU not found — skipping automatic include analysis")
    endif()
endif()

# ------------------------------------------------------------------------------
# SOURCES
# ------------------------------------------------------------------------------
set(astrotiger_src
   src/main.cpp
   src/fpe.cpp
   src/io.cpp
)

set(astrotiger_header
	include/autodiff.hpp
	include/cgs.hpp
	include/constants.hpp
	include/eos.hpp
	include/fpe.hpp
	include/gas_conserved.hpp
	include/gas_flux.hpp
	include/gas.hpp
	include/gas_primitive.hpp
	include/indices.hpp
	include/interval.hpp
	include/io.hpp
	include/left_right_state.hpp
	include/math.hpp
	include/matrix_fwd.hpp
	include/matrix.hpp
	include/multi_array.hpp
	include/opacity.hpp
	include/operators.hpp
	include/quantity.hpp
	include/rad_conserved.hpp
	include/rad_flux.hpp
	include/rational.hpp
	include/tensor.hpp
	include/units.hpp
	include/util.hpp
	include/vector.hpp
)

# ------------------------------------------------------------------------------
# CODEGEN: auto-build and run before astrotiger
# ------------------------------------------------------------------------------

# Assume submodule path is extern/codegen
add_subdirectory(extern/codegen)

# Location for the generated header
set(GENERATED_HEADER ${PROJECT_SOURCE_DIR}/include/dgTransforms.hpp)

# Define a custom command that runs codegen to create the header
add_custom_command(
    OUTPUT ${GENERATED_HEADER}
    COMMAND codegen codegen
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/include
    DEPENDS codegen
    COMMENT "Generating dgTransforms.hpp with codegen"
    VERBATIM
)

# Define a custom target representing the generated header
add_custom_target(generate_dgTransforms
    DEPENDS ${GENERATED_HEADER}
)

# Add the generated header to the sources so IDEs know about it
set(GENERATED_SOURCES ${GENERATED_HEADER})


file(GLOB_RECURSE CODEGEN_INPUTS
    ${PROJECT_SOURCE_DIR}/extern/codegen/src/*.cpp
    ${PROJECT_SOURCE_DIR}/extern/codegen/include/*.hpp
)

add_custom_command(
    OUTPUT ${GENERATED_HEADER}
    COMMAND codegen codegen
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/include
    DEPENDS codegen ${CODEGEN_INPUTS}
    COMMENT "Regenerating dgTransforms.hpp"
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

target_include_directories(astrotiger
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${GENERATED_SOURCE_DIR}
)

# Make astrotiger depend on this generated file
add_dependencies(astrotiger generate_dgTransforms)


set_property(TARGET astrotiger PROPERTY CUDA_ARCHITECTURES native)
set_property(TARGET astrotiger PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_compile_options(astrotiger PUBLIC
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O3 -g -fopt-info-vec-optimized >>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O0 -g -fno-omit-frame-pointer >>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -g -pg -O2 -DNDEBUG>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -G>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} --use_fast_math -DNDEBUG Xcompiler=-Ofast>>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -Xcompiler=-g -pg -DNDEBUG>>
)

target_link_options(astrotiger PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -rdynamic >>
)

# ------------------------------------------------------------------------------
# TESTS
# ------------------------------------------------------------------------------
add_subdirectory(extern/Catch2)

set(test_sources
  ${astrotiger_header}  
  src/fpe.cpp
  tests/test_autodiff.cpp
  tests/test_con2prim.cpp
  tests/test_implicit_radiation.cpp
  tests/test_matrix.cpp
  tests/test_radiation_eigensystem.cpp
  tests/test_srhd_eigensystem.cpp
)

add_executable(unit_tests ${test_sources})

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
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O0 -g -fsanitize=address,undefined -fno-omit-frame-pointer>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O3>>
)

target_link_options(unit_tests PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>: -fsanitize=address,undefined -fno-omit-frame-pointer>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:>>
)

include(CTest)
include(Catch)
catch_discover_tests(unit_tests)

