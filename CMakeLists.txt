cmake_minimum_required(VERSION 3.18)
project(astrotiger LANGUAGES CXX CUDA)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(ALL_FLAGS -Wall -Wno-comma-subscript -I${PROJECT_SOURCE_DIR}/include)
set(CXX_FLAGS  ${ALL_FLAGS} -fconstexpr-ops-limit=1000000000 -march=native)
set(CUDA_FLAGS ${ALL_FLAGS} -Xcompiler=-march=native)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(HPX REQUIRED)
set(astrotiger_src
   src/Main.cpp
   src/Options.cpp
)
set(astrotiger_header
   include/Approximate.hpp
   include/Box.hpp
   include/Child.hpp
   include/Definitions.hpp
   include/Face.hpp
   include/FixedPrecision.hpp
   include/Gas.hpp
   include/Integer.hpp
   include/IntegerPartition.hpp
   include/IO.hpp
   include/Leaf.hpp
   include/Math.hpp
   include/Matrix.hpp
   include/MortonKey.hpp
   include/Multidices.hpp
   include/Octree.hpp
   include/Options.hpp
   include/Point.hpp
   include/Rational.hpp
   include/Tensor.hpp
   include/YoungTableau.hpp
   include/Vector.hpp
)
add_hpx_executable(
  astrotiger
  COMPONENT_DEPENDENCIES
  SOURCES
    ${astrotiger_src}
  HEADERS
    ${astrotiger_header}
)
set_property(TARGET astrotiger PROPERTY CUDA_ARCHITECTURES native)
set_property(TARGET astrotiger PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_include_directories(astrotiger PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_compile_options(astrotiger PUBLIC
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -fPIE -DNDEBUG -Ofast>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -DNDEBUG -Xcompiler=-Ofast>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -fPIE -g>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -G -Xcompiler=-O0>>
)
target_link_options(astrotiger PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS}>>
)
