# PEnG-motif

 (C) Johannes Soeding, Markus Meier, Christian Roth

 [![Build Status](https://travis-ci.org/soedinglab/PEnG-motif.svg?branch=master)](https://travis-ci.org/soedinglab/PEnG-motif)
 [![License](https://img.shields.io/github/license/soedinglab/PEnG-motif.svg)](https://choosealicense.com/licenses/gpl-3.0/)
 [![Issues](https://img.shields.io/github/issues/soedinglab/PEnG-motif.svg)](https://github.com/soedinglab/PEnG-motif/issues)

PEnG-motif is an open-source software package for searching motifs (position specific weight matrices, PWMs) in a set of DNA sequences.

## Requirements

To compile and run PEnG-motif, you need
 * [CMake](http://cmake.org/) 2.8.12 or later
 * a recent C++ compiler (support for C++14)

### Optional requirements for running the scripts

  * python>=3.5
  * numpy
  * [BaMMmotif2](https://github.com/soedinglab/BaMMmotif2)


## Installation

### Cloning from GIT
If you want to compile the most recent version, simply clone the git repository.

	git clone https://github.com/soedinglab/PEnG-motif.git
	cd PEnG-motif

### Compilation
With the source code ready, simply run cmake with the default settings and libraries should be auto-detected:

	mkdir build
	cd build

Adjust ${HOME}/opt/PEnG if you want to change the installation directory

	cmake -DCMAKE_INSTALL_PREFIX=${HOME}/opt/PEnG ..
	make
	make install

### Environment setup
Add this line to your $HOME/.bashrc (or .zshrc...) to add peng_motif to your PATH:

	export PATH=${PATH}:${HOME}/opt/PEnG/bin

Update your environment:

	source $HOME/.bashrc


## Usage
PEnG-motif finds enriched sequences in a fasta file and writes the motifs to an output file in meme format.

The simplest invocation is:

```bash

  peng_motif <path/to/input.fasta> -o output.meme
```

For a list of all available options, please see the output of `peng_motif -h`.


## Interpreting the output

The PEnG-motif algorithm runs several phases; understanding the output printed to the console
can be important when interpreting the results.

### Phase 1: Counting base patterns
In the first phase, the occurrences of all 4^(pattern_length) kmers are counted. The expected number
of pattern counts are calculated under the assumption of a homogenous Markov model. From the observed
and expected counts a z-score is computed for each kmer.

```
pattern	       observed	     enrichment	         zscore

 CACTAG	            890	           2.40	          27.00
 GCCACC	           2456	           1.68	          26.08
 ```

 For each pattern that passes a predefined z-score threshold, the number of occurences, the enrichment of observed occurences over the expected and the calculated z-score are reported.

 ### Phase 2: Optimizing the base patterns
 In the second phase the selected base patterns from phase 1 are iteratively optimized by degenerating single nucleotides if the degenerated pattern achieves a higher score. By default we optimize a function based on [mutual information](https://en.wikipedia.org/wiki/Mutual_information) of observation and expectation.

 ```
 CACTAG	       890	 2.40	 -0.036173
 CWCTAG	      1687	 1.89	 -0.036744
optimization: CACTAG -> CWCTAG
```
For each iteration the current IUPAC pattern, the number of observed counts, the enrichment over the expected counts and the optimization score are printed.
Once the optimization runs into a local optimum, the original base pattern and its optimal IUPAC pattern are reported.

If the optimization runs into a pattern that has already been seen in a previous optimization, the optimization stops and the base pattern is removed.

At the end of the optimization phase, all IUPAC patterns, their occurrences, enrichment over the expected counts, and the calculated z-scores are reported.

```
pattern	       observed	     enrichment	         zscore

 CWCTAG	           1687	           1.89	          26.51
 GYCAYC	           4376	           1.55	          29.29
```

### Phase 3: selection, PWM generation

In phase 3, only the best scoring IUPAC PWMs are retained and are converted to PWMs.

### Phase 4: EM-optimization and merging
In the final phase, an expectation-maximization algorithm sharpens the PWMs. PWMs that have strong detectable overlaps are merged to form longer PWMs. The so optimized PWMs are written to the output file in [meme format](http://meme-suite.org/doc/meme-format.html).

## Tips and tricks

* PEnG-motif has many parameters that can be set to tune the results. The default parameters are suitable for most DNA binding factors. For RNA binding factors with short motifs a pattern length of 6 or 8 often yields better results.
* If your sequences are strand specific, pass the option `--strand PLUS` to avoid mixing the counts of reverse complemented base patterns.
* The memory requirement and runtime scales exponentially with the pattern length. It is *not* recommended to run with base pattern lengths larger than 12.

## Benchmarking
The `shoot_peng.py` python script can be used to annotate motifs according how well they can distinguish the given sequences from randomly generated ones.

Benchmarking requires [BaMMmotif2](https://github.com/soedinglab/BaMMmotif2).


## License

The PEnG-motif can be modified and distributed under the GPL-3.0 License.


## Acknowledgements

PEnG-motif uses shared code from BaMMmotif.
Many thanks to the developers of BaMMmotif!
