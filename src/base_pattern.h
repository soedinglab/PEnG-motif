/*
 * base_pattern.h
 *
 *  Created on: Jan 30, 2017
 *      Author: mmeier
 */

#ifndef SRC_BASE_PATTERN_H_
#define SRC_BASE_PATTERN_H_

#include <cstdlib>
#include <string>
#include "shared/SequenceSet.h"
#include "shared/Alphabet.h"
#include "shared/BackgroundModel.h"
#include "Global.h"


/**
    BasePattern deals with the encoding of Patterns
    The nucleotide encoding is defined in Alphabet

    The patterns are encoded in numerical values:
    (BaMM-)Alphabet: 'other'<->0, A<->1, C<->2, G<->3, T<->4
    (PEnG-)Alphabet: A<->0, C<->1, G<->2, T<->3
    Alphabet Size a = 4
    Pattern: "ATGC" <-> id: 0*a^0 + 3*a^1 + 2*a^2 + 1*a^3
*/
class BasePattern {
 public:
  //default de-/constructor
  BasePattern(const size_t pattern_length, Strand s, const int k, const int max_k,
              SequenceSet* sequence_set, BackgroundModel* bg);
  ~BasePattern();

  //inits factor and pattern_length
  void init(const size_t pattern_length);

  //returns factors
  size_t* getFactors();

  //returns pattern_length
  size_t getPatternLength();

  //get string to pattern id
  std::string toString(const size_t pattern_id);

  //get pattern id of reverse complementary pattern
  size_t getRevCompId(const size_t pattern_id);

  //get nucleotide id of pattern at pos
  int getNucleotideAtPos(const size_t pattern, const size_t pos);

  size_t getNumberPatterns();
  size_t* getPatternCounter();
  float* getBackgroundProb(const int order);
  size_t baseId2IUPACId(const size_t base_pattern);
  float getExpCountFraction(const size_t pattern, const size_t pseudo_expected_pattern_counts);
  float getLogPval(size_t pattern);
  float getOptimizationScore(const OPTIMIZATION_SCORE score_type, const size_t pattern, const size_t pseudo_expected_pattern_counts);
  size_t getLtot();

  void filter_base_patterns(const float zscore_threshold,
                                  const size_t count_threshold,
                                  std::vector<size_t>& selected_patterns);

 private:
  // position specific factors used for encoding of base patterns
  size_t* factor;

  // length of pattern
  size_t pattern_length;

  BackgroundModel* bg_model;
  size_t* pattern_counter;
  float** pattern_bg_probabilities;
  float* pattern_logp;
  float* pattern_zscore;

  size_t number_patterns;
  int max_k;
  int k;
  Strand strand;
  int alphabet_size;
  size_t n_sequences;
  size_t ltot;

  void count_patterns(SequenceSet* sequence_set);
  void count_patterns_minus_strand();
  size_t get_bg_id(const size_t pattern, const int curr_pattern_length, const int k);
  void calculate_bg_probabilities(BackgroundModel* model, const int alphabet_size, const int k, float* pattern_bg_probs);
  void calculate_bg_probability(float* background_model, const int alphabet_size,
                            const int k, int missing_pattern_length, size_t cur_pattern,
                            float cur_prob, float* final_probabilities);
  void calculate_log_pvalues(int ltot);
  void calculate_zscores(int ltot);
  float getMutualInformationScore(size_t pattern);
};

class sort_indices {
   private:
     float* mparr;
   public:
     sort_indices(float* parr) : mparr(parr) {}
     bool operator()(const size_t i, const size_t j) const { return mparr[i] > mparr[j]; }
};

#endif /* SRC_BASE_PATTERN_H_ */
