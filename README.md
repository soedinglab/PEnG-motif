# PEnG-motif

 (C) Johannes Soeding, Markus Meier, Christian Roth

 [![DOI](https://zenodo.org/badge/73262157.svg)](https://zenodo.org/badge/latestdoi/73262157)
 [![Build Status](https://travis-ci.org/soedinglab/PEnG-motif.svg?branch=master)](https://travis-ci.org/soedinglab/PEnG-motif)
 [![License](https://img.shields.io/github/license/soedinglab/PEnG-motif.svg)](https://choosealicense.com/licenses/gpl-3.0/)
 [![Issues](https://img.shields.io/github/issues/soedinglab/PEnG-motif.svg)](https://github.com/soedinglab/PEnG-motif/issues)


PEnG-motif is an open-source software package for searching motifs (position specific weight matrices, PWMs) in a set of DNA sequences.

As the core algorithm operates on kmers, the runtime is practically independent of the number and size of the input sequences. This makes PEnG-motif suitable for de-novo motif discovery on large sequence sets.

## Installation on macOS

The easiest way to get PEnG-motif on your mac is by downloading our precompiled binaries. If you choose that way, please download `peng_motif_macOS.zip` from our latest [release](https://github.com/soedinglab/PEnG-motif/releases) and unpack it.

After unzipping you will find a new binary `peng_motif` which you can call from command line by the full path to the binary `./peng_motif`, or directly by typing `peng_motif`, if you moved it to a location in your shell path.

## Installation on linux

The easiest way to get PEnG-motif on your 64bit linux computer is by downloading our precompiled binaries. If you choose that way, please download `peng_motif_linux_amd64.zip` from our latest [release](https://github.com/soedinglab/PEnG-motif/releases) and unpack it.

### Requirements for compiling from the scratch

 * [CMake](http://cmake.org/) 2.8.12 or later
 * gcc >=5.2 | gcc =4.9 (when manually setting CXXFLAG flag -std=c++1y)

### Installation procedure
Download the source code archive from our latest  [release](https://github.com/soedinglab/PEnG-motif/releases).

Unzip the source code, and navigate into the freshly unzipped folder in your terminal.

This code will compile and install PEnG-motif

```bash
INSTALL_DIR=/path/to/an/installation/directory/of/your/choice
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR
make && make install
```

You may want to add the path of `$INSTALL_DIR/bin` to your shell PATH, to simplify usage.

### Optional requirements for running the scripts
The repository also contains helper scripts - if you want to use them, you also need:
  * python>=3.5
  * numpy
  * [BaMMmotif2](https://github.com/soedinglab/BaMMmotif2)


## Using PEnG-motif
PEnG-motif finds enriched sequences in a fasta file and writes the motifs to an output file in meme format.

If the parent directory of of `peng_motif` is in your shell path, the simplest way to use PEng-motif is:

```bash

  peng_motif <path/to/input.fasta> -o output.meme
```

For a list of all available options, please see the output of `peng_motif -h`.

If you didn't put PEnG-motif in your shell path, you have to enter the full path to the peng binary, e.g. `/path/to/peng/peng_motif -h`.

## Interpreting the output

The PEnG-motif algorithm runs several phases; understanding the output printed to the console will help interpreting the final results.

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

## Docker image
PEnG-motif is included in our bamm-suite docker image `soedinglab/bamm-suite`.

After navigating into the directory with your fasta sequences use:

```
docker pull soedinglab/bamm-suite
docker run -v $UID:$GID -v $PWD/data soedinglab/bamm-suite peng_motif <fastafile.fa> <other options>

```

## License

The PEnG-motif can be modified and distributed under the GPL-3.0 License.


## Acknowledgements

PEnG-motif uses shared code from BaMMmotif.
Many thanks to the developers of BaMMmotif!
