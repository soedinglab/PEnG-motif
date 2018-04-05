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
  int bg_model_order = std::max(Global::bgModelOrder, Global::maxOptBgModelOrder);
  BackgroundModel* bgModel = new BackgroundModel(*Global::backgroundSequenceSet,
                          bg_model_order,
                          Global::bgModelAlpha,
                          Global::interpolateBG);

  #ifdef OPENMP
    omp_set_num_threads(Global::nr_threads);
  #endif

  //init peng with base patterns
  Peng peng(Global::strand, Global::bgModelOrder, Global::maxOptBgModelOrder,
                Global::inputSequenceSet, bgModel);


  //get merged degenerated iupac patterns from peng
  std::vector<IUPACPattern*> result;

  PengParameters params;
  params.max_pattern_length = Global::patternLength;
  params.zscore_threshold = Global::zscoreThreshold;
  params.count_threshold = Global::countThreshold;
  params.pseudo_counts = Global::pseudoCounts;
  params.opt_score_type = Global::optScoreType;
  // em options
  params.use_em = Global::useEm;
  params.em_saturation_factor = Global::emSaturationFactor;
  params.em_min_threshold = Global::emMinThreshold;
  params.em_max_iterations = Global::emMaxIterations;

  params.use_merging = Global::useMerging;
  params.bit_factor_merge_threshold = Global::mergeBitfactorThreshold;
  params.adv_pwm = Global::useAdvPWM;
  params.enrich_pseudocount_factor = Global::enrich_pseudocount_factor;
  params.minimum_processed_motifs = Global::minimum_processed_motifs;
  params.filter_neighbors = Global::filter_neighbors;
  params.max_optimized_patterns = Global::maximum_optimized_patterns;

  peng.process(params, result);

  peng.filter_redundancy(Global::mergeBitfactorThreshold, result);

  //print output
  if(Global::outputFilename) {
    peng.printShortMeme(result, Global::outputFilename, bgModel);
  }

  if(Global::jsonFilename) {
    peng.printJson(result, Global::jsonFilename, VERSION_NUMBER, bgModel);
  }

  //cleanup
  for(size_t i = 0; i < result.size(); i++) {
    delete result[i];
  }
  result.clear();

  delete bgModel;
}
