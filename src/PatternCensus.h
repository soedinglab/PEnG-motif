#ifndef PATTERN_CENSUS_H
#define PATTERN_CENSUS_H

#include "shared/SequenceSet.h"
#include "shared/Alphabet.h"
#include "shared/BackgroundModel.h"

class PatternCensus{
 public:
  PatternCensus(const int pattern_length, const int k, SequenceSet* sequence_set, BackgroundModel* bg);
  ~PatternCensus();

 private:
  void count_patterns(const int pattern_length, const int alphabet_size, SequenceSet* sequence_set);
  void count_patterns_reverse_strand(const int pattern_length, const int alphabet_size, size_t* pattern_counter);
  std::string getPatternFromNumber(size_t number);
  size_t get_bg_id(const size_t pattern, const int curr_pattern_length, const int k);
  void calculate_bg_probabilities(BackgroundModel* model, const int alphabet_size, const int k);
  void calculate_bg_probability(float* background_model, const int alphabet_size,
                            const int k, int missing_pattern_length, size_t pattern,
                            float cur_prob, float* final_probabilities);
  void calculate_log_pvalues(int ltot);
  void calculate_zscores(int ltot);
  void filter_nearest_neighbours(const int alphabet_size);

  size_t* pattern_counter;
  float* pattern_probabilities;
  float* pattern_logp;
  float* pattern_zscore;

  size_t number_patterns;
  int* factor;
  int pattern_length;
  int alphabet_size;
};

class sort_indices {
   private:
     size_t* mparr;
   public:
     sort_indices(size_t* parr) : mparr(parr) {}
     bool operator()(const size_t i, const size_t j) const { return mparr[i] > mparr[j]; }
};

#endif /* UTIL_H_ */
