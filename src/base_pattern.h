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
#include "iupac_pattern.h"
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

  friend class IUPACPattern;

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
  size_t getFastRevCompId(const size_t pattern_id);

  //get nucleotide id of pattern at pos
  int getNucleotideAtPos(const size_t pattern, const size_t pos);

  size_t getNumberPatterns();
  int getBackgroundOrder() const;
  float* getExpectedCounts() const;

  size_t* getPatternCounter();
  float* getBackgroundProb(const int order);
  float* getBackgroundProb();
  size_t baseId2IUPACId(const size_t base_pattern);
  float getExpCountFraction(const size_t pattern, const size_t pseudo_expected_pattern_counts);
  float getLogPval(size_t pattern);
  float getOptimizationScore(const OPTIMIZATION_SCORE score_type, const size_t pattern, const size_t pseudo_expected_pattern_counts);
  size_t getLtot();

  std::vector<size_t> select_base_patterns(const float zscore_threshold,
                                  const size_t count_threshold,
                                  bool single_stranded,
                                  bool filter_neighbors);

  std::vector<size_t> generate_double_stranded_em_optimization_patterns();
  std::unique_ptr<float[]> generate_correction_bg_probs(int order, BackgroundModel* model);
  auto generate_pattern_corrections(unsigned order, bool* seen_patterns, float* bg_probs);

  void print_patterns(std::vector<size_t> patterns);

  // helper methods
  inline size_t add_letter_to_the_right(size_t kmer, size_t position, int letter) {
    kmer += letter * factor[position];
    return kmer;
  }

  inline size_t get_bg_id(const size_t pattern, const int pattern_length, const int k) {

    /* this method extracts the rightmost (k+1)-mer from a kmer in PEnG representation
     * and returns its numerical value in BaMM representation
     *
     * PEnG: ACGT = 0*1 + 1*4 + 2*16 + 3*64
     * BaMM: ACGT = 0*64 + 1*16 + 2*4 + 3*1
    */
    size_t k_mer_pattern = 0;
    size_t* base_factors = BasePattern::getFactors();
    for(int i = pattern_length - k - 1; i < pattern_length; i++) {
      int c = BasePattern::getNucleotideAtPos(pattern, i);
      k_mer_pattern += c * base_factors[pattern_length - i - 1];
    }
    return k_mer_pattern;
  }

  inline size_t get_bg_id(const size_t pattern, const int pattern_length) {
    return get_bg_id(pattern, pattern_length, pattern_length - 1);
  }

 private:
  // position specific factors used for encoding of base patterns
  size_t* factor;

  // length of pattern
  size_t pattern_length;

  size_t* pattern_counter;
  float** pattern_bg_probabilities;
  float* pattern_logp;
  float* pattern_zscore;
  float* expected_counts;
  BackgroundModel* background_model;

  unsigned* half_revcomp;

  size_t number_patterns;
  int max_k;
  int k;
public:
  Strand getStrand() const;

private:
  Strand strand;
  int alphabet_size;
  size_t n_sequences;
  size_t ltot;

  void init_half_reverse_complements();
  void aggregate_double_strand_background();

  void count_patterns(SequenceSet* sequence_set);
  void count_patterns_single_strand(SequenceSet* sequence_set);
  void calculate_bg_probabilities(BackgroundModel* model, const int alphabet_size, const int k, float* pattern_bg_probs);
  void calculate_bg_probability(float* background_model, const int alphabet_size,
                            const int k, int missing_pattern_length, size_t cur_pattern,
                            float cur_prob, float* final_probabilities, size_t pattern_length);
  void calculate_log_pvalues();
  void calculate_zscores();
  void calculate_expected_counts();
  void calculate_expected_counts_single_stranded();
  float getMutualInformationScore(size_t pattern);
  void correct_counts();

};

class sort_indices {
   private:
     float* mparr;
   public:
     sort_indices(float* parr) : mparr(parr) {}
     bool operator()(const size_t i, const size_t j) const { return mparr[i] > mparr[j]; }
};

#endif /* SRC_BASE_PATTERN_H_ */
