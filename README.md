# PEnG-motif

 (C) Johannes Soeding, Markus Meier
 
The PEnG-motif an open-source software package for searching IUPAC patterns in a set of FASTA sequences.


## Requirements

To compile from source, you will need:
 * a recent C/C++ compiler (support for C++14)
 * [CMake](http://cmake.org/) 2.8.12 or later
 * For the shared code/scripts with [BaMMmotif](https://github.com/soedinglab/bamm-private)
 	* C++ library: [boost](http://www.boost.org/)
 	* R libraries from the CRAN package repository: fdrtool, zoo, argparse
 * From [BaMMmotif](https://github.com/soedinglab/bamm-private): plotAUSFC_rank.R and BaMMmotif need to be in the PATH


## Installation

### Cloning from GIT
If you want to compile the most recent version, simply clone the git repository. Then, from the repository root, initialize the BaMM submodule:

	git clone git@github.com:soedinglab/PEnG-motif.git
	cd PEnG-motif
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
For performing a search for overrepresented motives in a set of DNA sequences run PEnG-motif with the following command:

	scripts/shoot_peng.py -w 10 --pseudo-counts 10 -b 0.3 -a 1E3 -o output.meme input.fasta

You can get a detailed list of options for scripts/shoot_peng.py by calling without arguments.
This script is a wrapper for PEnG that will rerank the models using binaries and scripts from BaMM.


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

