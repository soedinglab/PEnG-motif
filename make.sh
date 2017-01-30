#!/usr/bin/zsh

rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX="/home/mmeier/opt/PEnG" ..
make -j 4
make install
cd ..
