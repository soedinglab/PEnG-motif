#ifndef PATTERN_CENSUS_H
#define PATTERN_CENSUS_H

#include "shared/SequenceSet.h"
#include "shared/Alphabet.h"
#include "shared/BackgroundModel.h"
#include <set>
#include <array>
#include <tuple>
#include "iupac_pattern.h"
#include "Global.h"

class Peng{
 public:
  Peng(const int pattern_length, Strand s, const int k, const int max_opt_k,
                SequenceSet* sequence_set, BackgroundModel* bg);
  ~Peng();
  void process(const float zscore_threshold, const size_t count_threshold, const int pseudo_counts,
                     const bool use_em, const float em_saturation_factor, const float min_em_threshold,
                     const int em_max_iterations, const bool use_merging, const float bit_factor_merge_threshold,
                     const bool adv_pwm,
                     std::vector<IUPACPattern*>& best_iupac_patterns);

  void printShortMeme(std::vector<IUPACPattern*>& best_iupac_patterns,
                      const std::string output_filename,
                      BackgroundModel* bg_model);

  void printJson(std::vector<IUPACPattern*>& best_iupac_patterns,
                 const std::string output_filename,
                 const std::string version_number,
                 BackgroundModel* bg_model);

 private:
  BackgroundModel* bg_model;
  size_t* pattern_counter;
  float** pattern_bg_probabilities;
  float* pattern_zero_bg_probabilities;
  float* pattern_logp;
  float* pattern_zscore;

  size_t number_patterns;
  int pattern_length;
  int max_k;
  int k;
  Strand strand;
  int alphabet_size;
  size_t ltot;

  void count_patterns(const int pattern_length, const int alphabet_size, const size_t number_patterns, SequenceSet* sequence_set, size_t* pattern_counter);

  void count_patterns_minus_strand(const int pattern_length, const int alphabet_size, const size_t number_patterns, size_t* pattern_counter);

  size_t get_bg_id(const size_t pattern, const int curr_pattern_length, const int k);

  void calculate_bg_probabilities(BackgroundModel* model, const int alphabet_size, const int k, float* pattern_bg_probs);

  void calculate_bg_probability(float* background_model, const int alphabet_size,
                            const int k, int missing_pattern_length, size_t cur_pattern,
                            float cur_prob, float* final_probabilities);

  void calculate_log_pvalues(int ltot);
  void calculate_zscores(int ltot);

  void filter_base_patterns(const int pattern_length, const int alphabet_size,
                                  const size_t number_patterns, const float zscore_threshold,
                                  const size_t count_threshold,
                                  float* pattern_zscore, std::vector<size_t>& selected_patterns);

  void optimize_iupac_patterns(std::vector<size_t>& selected_base_patterns,
                               std::vector<IUPACPattern*>& best_iupac_patterns);

  void filter_iupac_patterns(std::vector<IUPACPattern*>& iupac_patterns);

  void merge_iupac_patterns(const size_t pattern_length, const float bit_factor_merge_threshold,
                            BackgroundModel* bg, std::vector<IUPACPattern*>& iupac_patterns);

  void em_optimize_pwms(std::vector<IUPACPattern*>& iupac_patterns,
                                  const float saturation_factor,
                                  const float min_em_threshold, const int max_iterations,
                                  float* background_probabilities,
                                  std::vector<IUPACPattern*>& optimized_iupac_patterns);

  void calculate_prob_odds(const size_t pattern_length,
                      size_t curr_pattern, float curr_prob, int curr_length,
                      float** pwm, float* pattern_bg_probabilities, float* prob_odds);
};

class sort_indices {
   private:
     float* mparr;
   public:
     sort_indices(float* parr) : mparr(parr) {}
     bool operator()(const size_t i, const size_t j) const { return mparr[i] > mparr[j]; }
};

#endif /* UTIL_H_ */
