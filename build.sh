mkdir $1
cd $1
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON     \
      -DCMAKE_C_COMPILER_LAUNCHER=ccache     \
      -DCMAKE_CXX_COMPILER_LAUNCHER=ccache   \
      -DCMAKE_INSTALL_PREFIX=$HOME/local/$1/ \
      -DCMAKE_BUILD_TYPE=$1                  \
      ..
make -j


