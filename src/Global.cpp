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


char*	Global::alphabetType = nullptr;			              // alphabet type is defaulted to standard which is ACGT

char* Global::outputFilename = nullptr;                  // filename for IUPAC pattern output in short meme format
char* Global::jsonFilename = nullptr;                    // filename for IUPAC pattern output in json format

char* Global::inputSequenceFilename = nullptr;		        // filename with input FASTA sequences
char* Global::backgroundSequenceFilename = nullptr;      // filename with background FASTA sequences
SequenceSet* Global::inputSequenceSet = nullptr;         // input sequence Set
SequenceSet* Global::backgroundSequenceSet = nullptr;    // background sequence Set

OPTIMIZATION_SCORE Global::optScoreType = OPTIMIZATION_SCORE::MutualInformation;
float Global::enrich_pseudocount_factor = 0.005;

int Global::patternLength = 10;                        // length of patterns to be trained/searched
Strand Global::strand = Strand::BOTH_STRANDS;

bool Global::useEm = true;
float Global::emSaturationFactor = 1E4;
float Global::emMinThreshold = 0.08;
int Global::emMaxIterations = 10;

bool Global::useMerging = true;

int Global::pseudoCounts = 10;
bool Global::useAdvPWM = true;

float Global::zscoreThreshold = 10;
size_t Global::countThreshold = 1;
float Global::mergeBitfactorThreshold = 0.4;

// background model options
bool Global::interpolateBG = true;                    // calculate prior probabilities from lower-order probabilities
                                                      // instead of background frequencies of mononucleotides

int Global::bgModelOrder = 2;                         // background model order, defaults to 2
int Global::maxOptBgModelOrder = 2;
std::vector<float>  Global::bgModelAlpha( bgModelOrder+1, 1.0f );    // background model alpha

int Global::verbosity = 2;            	              // verbosity
int Global::nr_threads = 1;

bool Global::filter_neighbors = true;
unsigned Global::minimum_processed_motifs = 25;

void Global::init(int nargs, char* args[]){
  readArguments(nargs, args);

  Alphabet::init(alphabetType);

  // reverse complements are dealt with in peng directly. We read in single stranded sequences either way.
  inputSequenceSet = new SequenceSet(inputSequenceFilename, true);

  char* currBackgroundSequenceFilename;
  if(backgroundSequenceFilename) {
    currBackgroundSequenceFilename = backgroundSequenceFilename;
  }
  else {
    currBackgroundSequenceFilename = inputSequenceFilename;
  }

  backgroundSequenceSet = new SequenceSet(currBackgroundSequenceFilename, true);
}

void Global::readArguments(int nargs, char* args[]){
  if (nargs > 1 && !strcmp(args[1], "-h")) {
    Global::printHelp();
    exit(0);
  }
  else if (nargs > 1 && !strcmp(args[1], "-version")) {
    std::cout << "peng_motif version " << VERSION_NUMBER << std::endl;;
    exit(0);
  }

  if( nargs < 2 ) {			// At least input_file is required
    fprintf( stderr, "Error: Arguments are missing! \n" );
    printHelp();
    exit( -1 );
  }

  inputSequenceFilename = args[1];

  for (int i = 2; i < nargs; i++) {
    if (!strcmp(args[i], "-w")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following -w" << std::endl;
        exit(4);
      }
      patternLength = std::stoi(args[i]);
      if(patternLength % 2 == 1) {
        LOG(ERROR) << "Due to optimizations the pattern length has to be a multiple of 2" << std::endl;
        exit(4);
      }
    }
    else if (!strcmp(args[i], "--background-sequences")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following --background-sequences" << std::endl;
        exit(4);
      }
      backgroundSequenceFilename = args[i];
    }
    else if (!strcmp(args[i], "--optimization_score")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following --optimization_score" << std::endl;
        exit(4);
      }

      if(!strcmp(args[i], "LOGPVAL")) {
        optScoreType = OPTIMIZATION_SCORE::kLogPval;
      }
      else if(!strcmp(args[i], "ENRICHMENT")) {
        optScoreType = OPTIMIZATION_SCORE::kExpCounts;
      }
      else if(!strcmp(args[i], "MUTUAL_INFO")) {
        optScoreType = OPTIMIZATION_SCORE::MutualInformation;
      }
      else {
        printHelp();
        LOG(ERROR) << "Unknown expression following --optimization_score" << std::endl;
        exit(4);
      }
    }
    else if (!strcmp(args[i], "--enrich_pseudocount_factor")) {
        if (++i>=nargs) {
          printHelp();
          LOG(ERROR) << "No expression following --enrich_pseudocount_factor" << std::endl;
          exit(4);
        }
        enrich_pseudocount_factor = std::stof(args[i]);
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
    else if (!strcmp(args[i], "-j")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following -j" << std::endl;
        exit(4);
      }
      jsonFilename = new char[strlen(args[i]) + 1];
      strcpy(jsonFilename, args[i]);
    }
    else if (!strcmp(args[i], "-t")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following -t" << std::endl;
        exit(4);
      }
      zscoreThreshold = std::stof(args[i]);
    }
    else if (!strcmp(args[i], "--count-threshold")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following --count-threshold" << std::endl;
        exit(4);
      }
      countThreshold = std::stoi(args[i]);
    }
    else if (!strcmp(args[i], "-b")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following -b" << std::endl;
        exit(4);
      }
      mergeBitfactorThreshold = std::stof(args[i]);
    }
    else if (!strcmp(args[i], "--use-default-pwm")) {
      useAdvPWM = false;
    }
    else if (!strcmp(args[i], "--pseudo-counts")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following --pseudo-counts" << std::endl;
        exit(4);
      }
      pseudoCounts = std::stoi(args[i]);
    }
    else if (!strcmp(args[i], "--threads")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following -t" << std::endl;
        exit(4);
      }
      nr_threads = std::stoi(args[i]);
    }
    else if (!strcmp(args[i], "--no-em")) {
      useEm = false;
    }
    else if (!strcmp(args[i], "-a")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following -a" << std::endl;
        exit(4);
      }
      emSaturationFactor = std::stof(args[i]);
    }
    else if (!strcmp(args[i], "--em-threshold")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following --em-threshold" << std::endl;
        exit(4);
      }
      emMinThreshold = std::stof(args[i]);
    }
    else if (!strcmp(args[i], "--em-max-iterations")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following --em-max-iterations" << std::endl;
        exit(4);
      }
      emMaxIterations = std::stoi(args[i]);
    }
    else if (!strcmp(args[i], "--no-merging")) {
      useMerging = false;
    }
    else if (!strcmp(args[i], "--strand")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following --strand" << std::endl;
        exit(4);
      }

      if(!strcmp(args[i], "BOTH")) {
        strand = Strand::BOTH_STRANDS;
      }
      else if(!strcmp(args[i], "PLUS")) {
        strand = Strand::PLUS_STRAND;
      }
      else {
        printHelp();
        LOG(ERROR) << "Unknown expression following --strand" << std::endl;
        exit(4);
      }
    }
    else if (!strcmp(args[i], "--bg-model-order")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following --bg-model-order" << std::endl;
        exit(4);
      }
      bgModelOrder = std::stoi(args[i]);
    }
    else if (!strcmp(args[i], "--no-neighbor-filtering")) {
      filter_neighbors = false;
    }
    else if (!strcmp(args[i], "--minimum-processed-patterns")) {
      if (++i>=nargs) {
        printHelp();
        LOG(ERROR) << "No expression following --minimum-processed-patterns" << std::endl;
        exit(4);
      }
      minimum_processed_motifs = std::stoi(args[i]);
    }
    else if (!strcmp(args[i], "--version")) {
      std::cout << "peng_motif " << VERSION_NUMBER << std::endl;;
      exit(0);
    }
    else if (!strcmp(args[i], "-h")) {
      Global::printHelp();
      exit(0);
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
      "           best IUPAC motives will be written in OUTPUT_FILE\n"
      "           in minimal MEME format\n");
  printf("\n      -j, <OUTPUT_FILE>\n"
      "           best IUPAC motives will be written in OUTPUT_FILE\n"
      "           in JSON format\n");
  printf("\n      --background-sequences, <FASTA_FILE>\n"
      "           file with fasta sequences to be used for the\n"
      "           background model calculation\n");
  printf("\n      -t, <ZSCORE_THRESHOLD>\n"
      "           lower zscore threshold for basic patterns\n");
  printf("\n      -w, <PATTERN_LENGTH>\n"
      "           length of patterns to be searched\n");
  printf("\n      --bg-model-order, <BG_MODEL_ORDER>\n"
      "           order of the background model\n");
  printf("\n      --count-threshold, <COUNT_THRESHOLD>\n"
      "           lower threshold for counts of basic patterns\n");
  printf("\n      --strand, <PLUS|BOTH>\n"
      "           select the strands to work on\n");
  printf("\n      --optimization_score, <ENRICHMENT|LOGPVAL|MUTUAL_INFO>\n"
      "           select the iupac optimization score\n");
  printf("\n      --enrich_pseudocount_factor, <PSEUDO_COUNTS>\n"
	    "           add (enrich_pseudocount_factor x #seqs) pseudocounts\n"
	    "           in the EXPCOUNTS optimization\n");
  printf("\n      -b, <BIT_FACTOR_THRESHOLD>\n"
      "           bit factor threshold for merging IUPAC patterns\n");
  printf("\n      --no-em\n"
      "           shuts off the em optimization \n");
  printf("\n      -a, <EM_SATURATION_THRESHOLD>\n"
      "           saturation factor for em optimization \n");
  printf("\n      --em-threshold, <EM_THRESHOLD>\n"
      "           threshold for finishing the em optimization \n");
  printf("\n      --em-max-iterations, <EM_MAX_ITERATIONS>\n"
      "           max number of em optimization iterations\n");
  printf("\n      --no-merging\n"
      "           shuts off the merging \n");
  printf("\n      --use-default-pwm\n"
      "           use the default calculation of the pwm\n");
  printf("\n      --pseudo-counts, <PSEUDO_COUNTS>\n"
      "           number of pseudo-counts for optimization\n");
  printf("\n      --threads, <NUMBER_THREADS>\n"
      "           number of threads to be used for parallelization\n");
  printf("\n      --no-neighbor-filtering\n"
      "           do not filter similar base patterns before running the optimization\n");
  printf("\n      --minimum-processed-patterns <NUMBER_PATTERNS>\n"
      "           minimum number of iupac patterns that are selected for em optimization\n");
  printf("\n      --version\n"
      "           print the version number\n");
  printf("\n      -h\n"
      "           print this help \n");


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
