# PEnG-motif

 (C) Johannes Soeding, Markus Meier
 
The PEnG-motif an open-source software package for searching IUPAC patterns in a set of FASTA sequences.

## Requirements

To compile from source, you will need:
 * a recent C/C++ compiler
 * [CMake](http://cmake.org/) 2.8.12 or later

## Installation
We recommend compiling HHsuite on the machine that should run the computations so that it can be optimized for the appropriate CPU architecture.

### Packages
Some distributions incorporate HHsuite on their own:
* Ubuntu/Debian/etc. [DPKGs](http://packages.debian.org/source/sid/hhsuite) are provided by Laszlo Kajan.
* For Archlinux you can find a [PKGBUILD on aur](https://aur.archlinux.org/packages/hhsuite/)

### Release tarballs
The [release tarballs](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/releases/) should contain all required source files. Simply download and extract

### Cloning from GIT
If you want to compile the most recent version, simply clone the git repository. Then, from the repository root, initialize the ffindex submodule:

	git clone git@github.com:soedinglab/PEnG-motif.git
        cd PenG-motif
	git submodule init
	git submodule update


### Compilation
With the sourcecode ready, simply run cmake with the default settings and libraries should be auto-detected:

	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=${INSTALL_BASE_DIR} ..
	make
	make install


## Usage
For performing a search for common motives in a set of fasta sequences run PEnG-motif with the following command:

	peng_motif <input-fasta-file> -o <output-file>

You can get a detailed list of options for PEnG-motif by calling peng-motif without arguments.


## License

The PEnG-motif is distributed under the GPL-3.0 License.

## Notes

We are very grateful for bug reports! 
Please contact us at soeding@mpibpc.mpg.de

## Links

* [soeding lab](http://www.mpibpc.mpg.de/soeding)


## Acknowledgements
 
PEnG-motif uses shared code from BaMMmotif.
Many thanks to the developers of BaMMmotif for their great code!

