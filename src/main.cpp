/*
 * main.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: mmeier
 */

#include "log.h"
#include "Global.h"
#include "iupac_pattern.h"
#include "peng.h"

#ifdef OPENMP
  #include <omp.h>
#endif


int main(int nargs, char **args) {
  Global::init(nargs, args);

  //calculate background model
  BackgroundModel* bgModel = new BackgroundModel(*Global::backgroundSequenceSet,
                          Global::bgModelOrder,
                          Global::bgModelAlpha,
                          Global::interpolateBG);

  #ifdef OPENMP
    omp_set_num_threads(Global::nr_threads);
  #endif

  //init peng with base patterns
  Peng peng(Global::patternLength, Global::strand, Global::bgModelOrder,
                Global::inputSequenceSet, bgModel);


  //get merged degenerated iupac patterns from peng
  std::vector<IUPACPattern*> result;
  peng.process(Global::zscoreThreshold, Global::useEm,
               Global::emSaturationFactor, Global::emMinThreshold,
               Global::emMaxIterations, Global::useMerging, Global::mergeBitfactorThreshold,
               result);

  //print output
  if(Global::outputFilename != NULL) {
    peng.printShortMeme(result, Global::outputFilename, bgModel);
  }

  if(Global::jsonFilename != NULL) {
    peng.printJson(result, Global::jsonFilename, VERSION_NUMBER, bgModel);
  }

  //cleanup
  for(size_t i = 0; i < result.size(); i++) {
    delete result[i];
  }
  result.clear();

  delete bgModel;
}
