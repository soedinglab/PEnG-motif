/*
 * Global.cpp
 *
 *  Created on: Apr 4, 2016
 *      Author: wanwan
 */

#include <sys/stat.h>   			// get file status

#include "Global.h"
#include "log.h"
#include "shared/Alphabet.h"

#include <string.h>
#include <string>


char*	Global::alphabetType = NULL;			              // alphabet type is defaulted to standard which is ACGT

char* Global::outputFilename = NULL;                  // filename for IUPAC pattern output

char* Global::inputSequenceFilename = NULL;		        // filename with input FASTA sequences
SequenceSet* Global::inputSequenceSet = NULL;         // input sequence Set
SequenceSet* Global::backgroundSequenceSet = NULL;    // background sequence Set
bool Global::revcomp = false;                         // also search on reverse complement of sequences

int Global::patternLength = 8;                        // length of patterns to be trained/searched

float Global::zscoreThreshold = 1000;

// background model options
bool Global::interpolateBG = true;                    // calculate prior probabilities from lower-order probabilities
                                                      // instead of background frequencies of mononucleotides
int Global::bgModelOrder = 2;                         // background model order, defaults to 2
std::vector<float>  Global::bgModelAlpha( bgModelOrder+1, 1.0f );    // background model alpha

int Global::verbosity = 2;            	              // verbosity


void Global::init(int nargs, char* args[]){
	readArguments(nargs, args);

	Alphabet::init(alphabetType);

	inputSequenceSet = new SequenceSet(inputSequenceFilename);
	backgroundSequenceSet = new SequenceSet(inputSequenceFilename, true);
}

void Global::readArguments(int nargs, char* args[]){
	if( nargs < 2 ) {			// At least input_file is required
		fprintf( stderr, "Error: Arguments are missing! \n" );
		printHelp();
		exit( -1 );
	}

	inputSequenceFilename = args[1];

	for (int i = 2; i < nargs; i++) {
    if (!strcmp(args[i], "-W")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following -W" << std::endl;
        exit(4);
      }
      patternLength = std::stoi(args[i]);
    }
    else if (!strcmp(args[i], "-v")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following -v" << std::endl;
        exit(4);
      }
      verbosity = std::stoi(args[i]);
    }
    else if (!strcmp(args[i], "-o")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following -o" << std::endl;
        exit(4);
      }
      outputFilename = new char[strlen(args[i]) + 1];
      strcpy(outputFilename, args[i]);
    }
    else if (!strcmp(args[i], "-t")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following -t" << std::endl;
        exit(4);
      }
      zscoreThreshold = std::stof(args[i]);
    }
    else {
      LOG(WARNING) << "Ignoring unknown option " << args[i] << std::endl;
    }
	}

  alphabetType = new char[9];
  strcpy( alphabetType, "STANDARD");
}

void Global::printHelp(){
	printf("\n=================================================================\n");
	printf("\n Usage: peng_motif SEQFILE [options] \n\n");
	printf("\t SEQFILE: file with sequences in FASTA format. \n");
  printf("\n      -o, <OUTPUT_FILE>\n"
      "           best UIPAC motives will be written in OUTPUT_FILE.\n");
	printf("\n=================================================================\n");
}

void Global::destruct(){
	Alphabet::destruct();
	if( alphabetType )
	  delete[] alphabetType;
	if(outputFilename)
	  delete[] outputFilename;
	if(inputSequenceSet)
	  delete inputSequenceSet;
	if(backgroundSequenceSet)
	  delete backgroundSequenceSet;
	if(inputSequenceFilename)
	  delete[] inputSequenceFilename;
}
