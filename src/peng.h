#ifndef PATTERN_CENSUS_H
#define PATTERN_CENSUS_H

#include "shared/SequenceSet.h"
#include "shared/Alphabet.h"
#include "shared/BackgroundModel.h"
#include <set>
#include <array>
#include <vector>
#include <tuple>
#include "iupac_pattern.h"
#include "Global.h"

class PengParameters {
public:
	int max_pattern_length;
	float zscore_threshold;
	size_t count_threshold;
	int pseudo_counts;
	OPTIMIZATION_SCORE opt_score_type;
	float enrich_pseudocount_factor;

	// EM Parameters
	bool use_em;
	float em_saturation_factor;
	float em_min_threshold;
	int em_max_iterations;
	bool use_merging;
	float bit_factor_merge_threshold;
	bool adv_pwm;
	unsigned minimum_processed_motifs;
  bool filter_neighbors;
  int max_optimized_patterns;
};

class Peng{
 public:
  Peng(Strand s, const int k, const int max_opt_k,
                SequenceSet* sequence_set, BackgroundModel* bg);
  ~Peng();
  void process(PengParameters& params, std::vector<IUPACPattern*>& best_iupac_patterns);

  void filter_redundancy(const float merge_bit_factor_threshold, std::vector<IUPACPattern*>& iupac_patterns);

  void printShortMeme(std::vector<IUPACPattern*>& best_iupac_patterns,
                      const std::string output_filename,
                      BackgroundModel* bg_model);

  void printJson(std::vector<IUPACPattern*>& best_iupac_patterns,
                 const std::string output_filename,
                 const std::string version_number,
                 BackgroundModel* bg_model);

 private:
  BackgroundModel* bg_model;
  SequenceSet* sequence_set;
  int max_k;
  int k;
  int alphabet_size;
  size_t n_sequences;
  Strand strand;


  void optimize_iupac_patterns(OPTIMIZATION_SCORE score_type,
                               BasePattern* base_patterns,
                               std::vector<size_t>& selected_base_patterns,
                               std::vector<IUPACPattern*>& best_iupac_patterns,
							   float enrich_pseudocount_factor);

  void filter_iupac_patterns(const size_t pattern_length, unsigned minimum_retained_motifs, std::vector<IUPACPattern*>& iupac_patterns);

  void merge_iupac_patterns(const size_t pattern_length, const float bit_factor_merge_threshold,
                            BackgroundModel* bg, std::vector<IUPACPattern*>& iupac_patterns);

  void em_optimize_pwms(std::vector<IUPACPattern*>& iupac_patterns,
                        BasePattern* base_patterns,
                                  const float saturation_factor,
                                  const float min_em_threshold, const int max_iterations,
                                  float* background_probabilities,
                                  std::vector<IUPACPattern*>& optimized_iupac_patterns);

  void calculate_prob_odds(const size_t pattern_length,
                      size_t curr_pattern, float curr_prob, int curr_length,
                      float** pwm, float* pattern_bg_probabilities, size_t* factors, float* prob_odds);
};

#endif /* UTIL_H_ */
