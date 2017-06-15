#!/bin/zsh

# change compilation to clang:
# http://google-engtools.blogspot.de/2011/05/c-at-google-here-be-dragons.html
export CC=clang
export CXX=clang++

rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX="$HOME/opt/PEnG" ..

# run automatic bugfinder for clang: http://clang-analyzer.llvm.org/
# -B forces to build again... needed for the bug tracking
scan-build make -B

make -j 4
make install
cd ..
