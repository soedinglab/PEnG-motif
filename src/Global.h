/*
 * Global.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include "shared/SequenceSet.h"

class Global{

public:
	static char* alphabetType;						       // provide alphabet type
	static char* outputFilename;                 // filename for IUPAC pattern output
	static char* inputSequenceFilename;				   // filename of positive sequence FASTA file
	static SequenceSet*	inputSequenceSet;				 // positive Sequence Set
	static SequenceSet* backgroundSequenceSet;   // background Sequence Set
	static bool revcomp;                         // also search on reverse complement of sequences

	static int patternLength;                    // length of pattern to be searched/trained

	static int bgModelOrder;						         // background model order, defaults to 2
	static bool interpolateBG;                   // calculate prior probabilities from lower-order probabilities
	                                             // instead of background frequencies of mononucleotides
	static std::vector<float> bgModelAlpha;      // background model alpha


	static int verbosity;							           // verbose printouts, defaults to false

	static void init( int nargs, char* args[] );
	static void destruct();

private:
	static void readArguments( int nargs, char* args[] );
	static void printHelp();
};


#endif /* GLOBAL_H_ */
