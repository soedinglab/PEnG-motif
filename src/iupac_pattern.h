/*
 * iupac_pattern.h
 *
 *  Created on: Jan 30, 2017
 *      Author: mmeier
 */

#ifndef SRC_IUPAC_PATTERN_H_
#define SRC_IUPAC_PATTERN_H_

#include <cstdlib>
#include <string>
#include <vector>
#include "Global.h"


class BasePattern;

static const int MIN_MERGE_OVERLAP = 6;


/**
    IUPACPattern captures degenerated (IUPAC) patterns
    The nucleotides are encoded in iupac_alphabet

    The iupac pattern encoding is a "binary" representation:
    Alphabet: A<->0, C<->1, G<->2, T<->3, S<->4, W<->5, R<->6, Y<->7, M<->8, K<->9, N<->10
    Alphabet size, a = 11
    Pattern: ATGCNW <-> 0*a^0 + 3*a^1 + 2*a^2 + 1*a^3 + 10*a^4 + 6*a^5

    An less limited representation of degenerated patterns is the pwm
*/
class IUPACPattern {
 public:
  static size_t* iupac_factor;

  static void init(size_t pattern_length, float* bg_model);
  static std::string toString(size_t pattern_id, size_t pattern_length);
  static int getNucleotideAtPos(const size_t pattern, const size_t pos);
  static size_t baseToId(const size_t base_pattern, const size_t pattern_length);
  static size_t toId(std::string base_pattern, const size_t pattern_length);
  static std::tuple<float, int, bool> calculate_S(IUPACPattern* p1, IUPACPattern* p2, Strand s, float* background);
  static void normalize_pwm(const int pattern_length, float** pwm);
  static float calculate_s(float** p1_pwm, float** p2_pwm, float* background, const int offset1, const int offset2, const int l);

  static void initIUPACProfile(const float mixing_factor, const float mixin_bias, float* bg_model);

  IUPACPattern(size_t iupac_pattern, size_t pattern_length);
  IUPACPattern(IUPACPattern* ori, float** pwm);
  IUPACPattern(IUPACPattern* longer_pattern, IUPACPattern* shorter_pattern,
               bool is_comp, float* background, const int shift);
  ~IUPACPattern();

  size_t get_pattern();
  std::string get_pattern_string();
  int get_pattern_length();

  float getExpCountFraction(const size_t pseudo_expected_pattern_counts);
  float getLogPval();
  float getMutualInformationScore(unsigned int n_sequences);
  float getOptimizationScore(OPTIMIZATION_SCORE score_type, const size_t pseudo_expected_pattern_counts,
                             unsigned int n_sequences);

  float get_bg_p();

  float** get_pwm();
  float** get_comp_pwm();

  void calculate_comp_pwm();
  void update_pwm(float** new_pwm);

  size_t get_sites();
  size_t* get_local_sites();

  float getExpectedCounts() const;
  float getZscore() const;

  std::vector<size_t>& get_base_patterns();
  int get_optimization_bg_model_order();
  void set_optimization_bg_model_order(int order);

  void calculate_pwm(BasePattern* base_pattern, const int pseudo_counts, size_t* pattern_counter, float* background_model);
  void calculate_adv_pwm(BasePattern* base_pattern, const int pseudo_counts, size_t* pattern_counter, float* background_model);
  void aggregate_attributes_from_basepatterns(BasePattern*);


  bool operator<(const IUPACPattern& rhs) const;

  std::vector<size_t> generate_base_patterns(BasePattern* basepatterns, size_t iupac_pattern);

  unsigned long count_combined_occurences(BasePattern*, size_t iupac_pattern);
  std::vector<size_t> basepatterns_from_iupac_single_stranded(BasePattern*, size_t iupac_pattern);
  std::vector<size_t> basepatterns_from_iupac_double_stranded(BasePattern*, size_t iupac_pattern);

  static void find_base_patterns(BasePattern* base_pattern, const size_t pattern, const size_t pattern_length, std::vector<size_t>& base_patterns);

 private:
  static float calculate_d(float** p1_pwm, float** p2_pwm, const int offset1, const int offset2, const int l, const float epsilon = 1E-4);
  static float calculate_d_bg(float** p_pwm, float* background, const int l, const int offset = 0, const float epsilon = 1E-4);

  float calculate_merged_pvalue(IUPACPattern* longer_pattern, IUPACPattern* shorter_pattern, bool is_comp, float* background, const int shift);

  static float** iupac_profile;
  static float* log_bonferroni;
  size_t pattern_length;

  size_t pattern;
  float log_pvalue;
  float zscore;
  float bg_p;
  float expected_counts;
  int optimization_bg_model_order;
  size_t n_sites;
  size_t* local_n_sites;
  float** pwm;
  float** comp_pwm;
  std::vector<size_t> base_patterns;
  bool merged;
};

bool sort_IUPAC_patterns(IUPACPattern* a, IUPACPattern* b);



#endif /* SRC_IUPAC_PATTERN_H_ */
