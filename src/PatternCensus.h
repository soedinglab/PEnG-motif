#ifndef PATTERN_CENSUS_H
#define PATTERN_CENSUS_H

#include "shared/SequenceSet.h"
#include "shared/Alphabet.h"
#include "shared/BackgroundModel.h"
#include <set>

class PatternCensus{
 public:
  PatternCensus(const int pattern_length, const int k, SequenceSet* sequence_set, BackgroundModel* bg);
  ~PatternCensus();

 private:

  void count_patterns(const int pattern_length, const int alphabet_size, SequenceSet* sequence_set);
  void count_patterns_minus_strand(const int pattern_length, const int alphabet_size, size_t* pattern_counter);
  std::string getPatternFromNumber(size_t number);
  std::string getIUPACPatternFromNumber(size_t pattern_id);
  size_t get_bg_id(const size_t pattern, const int curr_pattern_length, const int k);
  void calculate_bg_probabilities(BackgroundModel* model, const int alphabet_size, const int k);
  void calculate_bg_probability(float* background_model, const int alphabet_size,
                            const int k, int missing_pattern_length, size_t pattern,
                            float cur_prob, float* final_probabilities);
  void calculate_log_pvalues(int ltot);
  void calculate_zscores(int ltot);
  void filter_nearest_neighbours(const int alphabet_size, std::set<size_t>& selected_patterns);
  size_t get_rev_pattern_id(const size_t pattern_id, const int pattern_length,
                            const int* factor, const int* rev_factor);
  float find_base_patterns(const size_t mutated_pattern, const int pattern_length,
                           int* factor, int* iupac_factor, std::set<size_t>& base_patterns);

  float calculate_logpvalue_of_iupac_pattern(size_t mutated_pattern, const int ltot,
                                            float* base_background_prob, size_t* base_counts);

  void optimize_iupac_patterns(const int pattern_length, std::set<size_t>& selected_patterns,
                               std::set<size_t>& best_iupac_patterns);

  void filter_iupac_patterns(std::set<size_t>& iupac_patterns,
                             const int pattern_length, int* iupac_factor);

  void get_pwm_for_iupac_patterns(std::set<size_t>& best_iupac_patterns,
                                  const int pattern_length, size_t* pattern_counter);



  size_t* pattern_counter;
  float* pattern_bg_probabilities;
  float* pattern_logp;
  float* pattern_zscore;

  size_t number_patterns;
  int* factor;
  int* rev_factor;
  int* iupac_factor;
  float* log_bonferroni;
  int pattern_length;
  int alphabet_size;
  int ltot;
};

class sort_indices {
   private:
     float* mparr;
   public:
     sort_indices(float* parr) : mparr(parr) {}
     bool operator()(const size_t i, const size_t j) const { return mparr[i] > mparr[j]; }
};

#endif /* UTIL_H_ */
