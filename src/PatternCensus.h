#ifndef PATTERN_CENSUS_H
#define PATTERN_CENSUS_H

#include "shared/SequenceSet.h"
#include "shared/Alphabet.h"
#include "shared/BackgroundModel.h"
#include <set>
#include <array>
#include <tuple>

class BasePattern{
 public:
  static void init(size_t pattern_length);
  static size_t* factor;
  static size_t* rev_factor;
  static size_t pattern_length;

  static std::string getPatternFromNumber(size_t pattern_id);
  static size_t get_rev_pattern_id(const size_t pattern_id);

  static int getNucleotideAtPos(const size_t pattern, const size_t pos);
};

class IUPACPattern {
 public:
  static size_t* iupac_factor;

  static void init(size_t pattern_length);
  static std::string getIUPACPatternFromNumber(size_t pattern_id, size_t pattern_length);
  static int getNucleotideAtPos(const size_t pattern, const size_t pos);
  static size_t getIUPACPattern(const size_t base_pattern, const size_t pattern_length);
  static size_t getIUPACPattern(std::string base_pattern, const size_t pattern_length);

  static float calculate_s(IUPACPattern& p1, IUPACPattern& p2, float* background, const int offset1, const int offset2, const int l);
  static float calculate_d(IUPACPattern& p1, IUPACPattern& p2, const int offset1, const int offset2, const int l);
  static float calculate_d_bg(IUPACPattern& p, float* background, const int l, const int offset = 0);
  static std::tuple<float, int> calculate_S(IUPACPattern* p1, IUPACPattern* p2, float* background);

  float calculate_merged_pvalue(IUPACPattern* longer_pattern, IUPACPattern* shorter_pattern, float* background, const int shift);

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

  void find_base_patterns(const size_t pattern, const size_t pattern_length, std::vector<size_t>& base_patterns);
  void count_sites(size_t* pattern_counter);
  void calculate_pwm(size_t* pattern_counter);
  void calculate_adv_pwm(size_t* pattern_counter);
  void calculate_log_pvalue(const int ltot,
                                            float* base_background_prob,
                                            size_t* base_counts);

  bool operator<(const IUPACPattern& rhs) const;

 private:
  static float* log_bonferroni;
  size_t pattern_length;

  size_t pattern;
  float log_pvalue;
  float bg_p;
  size_t n_sites;
  size_t* local_n_sites;
  float** pwm;
  std::vector<size_t> base_patterns;
};

class PatternCensus{
 public:
  PatternCensus(const int pattern_length, const int k, const float zscore_threshold,
                SequenceSet* sequence_set, BackgroundModel* bg, const char* outputFilename, const char* jsonFilename);
  ~PatternCensus();

 private:
  void count_patterns(const int pattern_length, const int alphabet_size, SequenceSet* sequence_set);
  void count_patterns_minus_strand(const int pattern_length, const int alphabet_size, size_t* pattern_counter);

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

  void optimize_iupac_patterns(const float zscore_threshold,
                               std::set<size_t>& selected_base_patterns,
                               std::vector<IUPACPattern*>& best_iupac_patterns);

  void filter_iupac_patterns(std::vector<IUPACPattern*>& iupac_patterns);

  void printShortMeme(std::vector<IUPACPattern*>& best_iupac_patterns,
                      const std::string output_filename,
                      const std::string version_number,
                      BackgroundModel* bg_model);

  void printJson(std::vector<IUPACPattern*>& best_iupac_patterns,
                 const std::string output_filename,
                 const std::string version_number,
                 BackgroundModel* bg_model);


  size_t* pattern_counter;
  float* pattern_bg_probabilities;
  float* pattern_logp;
  float* pattern_zscore;

  size_t number_patterns;
  int pattern_length;
  int alphabet_size;
  size_t ltot;
};

bool sort_IUPAC_patterns(IUPACPattern* a, IUPACPattern* b);

class sort_indices {
   private:
     float* mparr;
   public:
     sort_indices(float* parr) : mparr(parr) {}
     bool operator()(const size_t i, const size_t j) const { return mparr[i] > mparr[j]; }
};

#endif /* UTIL_H_ */
