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

  static void init(size_t pattern_length);
  static std::string toString(size_t pattern_id, size_t pattern_length);
  static int getNucleotideAtPos(const size_t pattern, const size_t pos);
  static size_t baseToId(const size_t base_pattern, const size_t pattern_length);
  static size_t toId(std::string base_pattern, const size_t pattern_length);
  static std::tuple<float, int> calculate_S(IUPACPattern* p1, IUPACPattern* p2, float* background);

  IUPACPattern(size_t iupac_pattern, size_t pattern_length);
  IUPACPattern(IUPACPattern* longer_pattern, IUPACPattern* shorter_pattern, float* background, const int shift);
  ~IUPACPattern();

  size_t get_pattern();
  int get_pattern_length();
  float get_log_pvalue();
  float get_bg_p();
  float** get_pwm();
  size_t get_sites();
  size_t* get_local_sites();
  std::vector<size_t>& get_base_patterns();

  void count_sites(size_t* pattern_counter);
  void calculate_pwm(size_t* pattern_counter);
  void calculate_adv_pwm(size_t* pattern_counter, float* background_model);
  void calculate_log_pvalue(const int ltot,
                            float* base_background_prob,
                            size_t* base_counts);

  bool operator<(const IUPACPattern& rhs) const;

 private:
  static float calculate_s(IUPACPattern& p1, IUPACPattern& p2, float* background, const int offset1, const int offset2, const int l);
  static float calculate_d(IUPACPattern& p1, IUPACPattern& p2, const int offset1, const int offset2, const int l);
  static float calculate_d_bg(IUPACPattern& p, float* background, const int l, const int offset = 0);
  static void find_base_patterns(const size_t pattern, const size_t pattern_length, std::vector<size_t>& base_patterns);

  float calculate_merged_pvalue(IUPACPattern* longer_pattern, IUPACPattern* shorter_pattern, float* background, const int shift);

  static float* log_bonferroni;
  size_t pattern_length;

  size_t pattern;
  float log_pvalue;
  float bg_p;
  size_t n_sites;
  size_t* local_n_sites;
  float** pwm;
  std::vector<size_t> base_patterns;
  bool merged;
};

bool sort_IUPAC_patterns(IUPACPattern* a, IUPACPattern* b);



#endif /* SRC_IUPAC_PATTERN_H_ */
