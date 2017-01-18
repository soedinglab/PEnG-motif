/*
 * main.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: mmeier
 */

#include "log.h"
#include "Global.h"
#include "PatternCensus.h"

#ifdef OPENMP
  #include <omp.h>
#endif


int main(int nargs, char **args) {
  Global::init(nargs, args);

  BackgroundModel* bgModel = new BackgroundModel(*Global::backgroundSequenceSet,
                          Global::bgModelOrder,
                          Global::bgModelAlpha,
                          Global::interpolateBG );

  #ifdef OPENMP
    omp_set_num_threads(Global::nr_threads);
  #endif

  PatternCensus(Global::patternLength, Global::bgModelOrder,
                Global::zscoreThreshold, Global::inputSequenceSet,
                bgModel, Global::outputFilename, Global::jsonFilename);

  delete bgModel;
}
