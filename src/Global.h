/*
 * Global.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include "shared/SequenceSet.h"

enum class Strand {
  PLUS_STRAND, BOTH_STRANDS
};

enum class OPTIMIZATION_SCORE {
	kLogPval = 0,
	kExpCounts = 1,
	MutualInformation
};

const std::string VERSION_NUMBER("1.0.0");

class Global{
public:
	static char* alphabetType;						       // provide alphabet type
	static char* outputFilename;                 // filename for IUPAC pattern output in short meme format
	static char* jsonFilename;                   // filename for IUPAC pattern output in json format
	static char* inputSequenceFilename;				   // filename of positive sequence FASTA file
	static char* backgroundSequenceFilename;     // filename of background sequence FASTA file
	static SequenceSet*	inputSequenceSet;				 // positive Sequence Set
	static SequenceSet* backgroundSequenceSet;   // background Sequence Set
	static bool revcomp;                         // also search on reverse complement of sequences

	static OPTIMIZATION_SCORE optScoreType;
	static float enrich_pseudocount_factor;

	static int patternLength;                    // length of pattern to be searched/trained
	static float zscoreThreshold;
	static size_t countThreshold;
	static Strand strand;

	static bool useEm;
	static float emSaturationFactor;
	static float emMinThreshold;
	static int emMaxIterations;

	static bool useMerging;

	static float mergeBitfactorThreshold;

	static bool useAdvPWM;
	static int pseudoCounts;

	static int bgModelOrder;						         // background model order, defaults to 2
	static int maxOptBgModelOrder;               // max background model order for optimization, defaults to 3
	static bool interpolateBG;                   // calculate prior probabilities from lower-order probabilities
	                                             // instead of background frequencies of mononucleotides
	static std::vector<float> bgModelAlpha;      // background model alpha

	static int nr_threads;

	static int verbosity;							           // verbose printouts, defaults to false

	static void init( int nargs, char* args[] );
	static void destruct();

	static bool filter_neighbors;
	static unsigned minimum_processed_motifs;

private:
	static void readArguments( int nargs, char* args[] );
	static void printHelp();
};


#endif /* GLOBAL_H_ */
