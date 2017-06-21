# PEnG-motif

 (C) Johannes Soeding, Markus Meier

PEnG-motif is an open-source software package for searching motifs (position specific weight matrices, PWMs) in a set of DNA sequences.
It is intended to be a fast prefilter for a future release of [BaMMmotif](https://github.com/soedinglab/BaMMmotif).


## Requirements

To compile from source, you will need:
 * [CMake](http://cmake.org/) 2.8.12 or later
 * C++ compiler (support for C++14)


## Installation

### Cloning from GIT
If you want to compile the most recent version, simply clone the git repository.

	git clone git@github.com:soedinglab/PEnG-motif.git
	cd PEnG-motif

### Compilation
With the source code ready, simply run cmake with the default settings and libraries should be auto-detected:

	mkdir build
	cd build

Adjust ${HOME}/opt/PEnG if you want to change the installation directory

	cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=${HOME}/opt/PEnG ..
	make
	make install

### Environment setup
Add this line to your $HOME/.bashrc (or .zshrc...) to add peng_motif to your PATH:

	export PATH=${PATH}:${HOME}/opt/PEnG/bin

Update your environment:

	source $HOME/.bashrc


## Usage
For performing a search for overrepresented motives in a set of DNA sequences run PEnG-motif with the following command:

  peng_motif -w 10 --pseudo-counts 10 -b 0.4 -a 1E4 -o output.meme input.fasta

You can get a detailed list of options for peng_motif by calling it without arguments.
This script is a wrapper for PEnG that will rerank the models using binaries and scripts from BaMM.

## Future
With a future release of [BaMMmotif](https://github.com/soedinglab/BaMMmotif) we can rerank
the models of PEnG with a proper score. For this purpose the python script shoot_peng.py
exists.


## License

The PEnG-motif is distributed under the GPL-3.0 License.


## Notes

We are very grateful for bug reports!
Please contact us at soeding@mpibpc.mpg.de


## Links

* [soeding lab](http://www.mpibpc.mpg.de/soeding)


## Acknowledgements

PEnG-motif uses shared code from BaMMmotif.
Many thanks to the developers of BaMMmotif!
