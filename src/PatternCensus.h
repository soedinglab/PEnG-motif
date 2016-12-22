#ifndef PATTERN_CENSUS_H
#define PATTERN_CENSUS_H

#include "shared/SequenceSet.h"
#include "shared/Alphabet.h"
#include "shared/BackgroundModel.h"
#include <set>
#include <array>

class BasePattern{
 public:
  static void init(size_t pattern_length);
  static int* factor;
  static int* rev_factor;
  static size_t pattern_length;

  static std::string getPatternFromNumber(size_t pattern_id);
  static size_t get_rev_pattern_id(const size_t pattern_id);

  static int getNucleotideAtPos(const size_t pattern, const size_t pos);
};

class IUPACPattern {
 public:
  static int* iupac_factor;

  static void init(size_t pattern_length);
  static std::string getIUPACPatternFromNumber(size_t pattern_id);
  static int getNucleotideAtPos(const size_t pattern, const size_t pos);
  static size_t getIUPACPattern(size_t base_pattern);
  static size_t getIUPACPattern(std::string base_pattern);

  IUPACPattern(size_t iupac_pattern);
  ~IUPACPattern();

  size_t get_pattern();
  float get_log_pvalue();
  float get_bg_p();
  float** get_pwm();
  size_t get_sites();

  void find_base_patterns();
  void count_sites(size_t* pattern_counter);
  void calculate_pwm(size_t* pattern_counter);
  void calculate_logpvalue_of_iupac_pattern(const int ltot,
                                            float* base_background_prob,
                                            size_t* base_counts);

  bool operator<(const IUPACPattern& rhs) const;

 private:
  static float* log_bonferroni;
  static size_t pattern_length;

  size_t pattern;
  float log_pvalue;
  float bg_p;
  size_t n_sites;
  float** pwm;
  std::set<size_t> base_patterns;
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
                               std::set<IUPACPattern*>& best_iupac_patterns);

  void filter_iupac_patterns(std::set<IUPACPattern*>& iupac_patterns);

  void printShortMeme(std::set<IUPACPattern*>& best_iupac_patterns,
                      const std::string output_filename,
                      const std::string version_number,
                      BackgroundModel* bg_model);

  void printJson(std::set<IUPACPattern*>& best_iupac_patterns,
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
  int ltot;
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
