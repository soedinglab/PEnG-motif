# macOS version requires homebrew's gcc@7 and cmake

OS_NAME=$(uname -s)

if [ "${OS_NAME}" == "Darwin" ]; then 
	export CXX=g++-7
	export CC=gcc-7
	export LDFLAGS="-static-libgcc -static-libstdc++"
elif [ "${OS_NAME}" == "Linux" ]; then 
	export LDFLAGS="-static -static-libgcc -static-libstdc++"
fi

BUILD_DIR=build_release
rm -rf $BUILD_DIR
mkdir $BUILD_DIR

cd $BUILD_DIR
cmake ..
make -j4

cd bin

if [ "${OS_NAME}" == "Darwin" ]; then
	zip peng_motif_macOS.zip peng_motif
elif [ "${OS_NAME}" == "Linux" ]; then
	zip peng_motif_linux_amd64.zip peng_motif
fi

mv *.zip ../
