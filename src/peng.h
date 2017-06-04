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
#include "gap_mask.h"

class Peng{
 public:
  Peng(const int pattern_length, Strand s, const int k,
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
  std::vector<GapMask*> masks;

  SequenceSet* sequence_set;
  int k;

  BackgroundModel* bg_model;
  size_t* pattern_counter;
  float* pattern_bg_probabilities;
  float* pattern_logp;
  float* pattern_zscore;

  size_t number_patterns;
  int pattern_length;
  Strand strand;
  int alphabet_size;
  size_t ltot;

  void count_patterns(const int pattern_length, const int alphabet_size, const size_t number_patterns, GapMask* mask, SequenceSet* sequence_set, size_t* pattern_counter);

  /**
      Helper function for the background model
      Extracts a k_mer id from a pattern to

      @param pattern is a base pattern id
      @param curr_pattern_length is the current position in the pattern
      @param k is the length of the kmer
      @return the id of the k_mer ending at curr_pattern_pos
  */
  void count_patterns_minus_strand(const int pattern_length, const int alphabet_size, const size_t number_patterns, size_t* pattern_counter);

  /**
      Helper function for the background model
      Extracts a k_mer id from a pattern to

      @param pattern is a base pattern id
      @param curr_pattern_length is the current position in the pattern
      @param k is the length of the kmer
      @return the id of the k_mer ending at curr_pattern_pos
  */
  size_t get_bg_id(const size_t pattern, const int curr_pattern_length, const int k, size_t* factors);



  /**
      Calculates the background probabilities of base patterns

      @param model the background model
      @param alphabet_size the size of the alphabet
      @param k length of mers to be used in the background model
  */
  void calculate_bg_probabilities(const int alphabet_size, const int k, GapMask* mask, BackgroundModel* model, size_t* factors, float* full_mask_patterns);

  /**
      Recursive calculation of background probabilities of base patterns

      @param background_model the background model
      @param alphabet_size the size of the alphabet
      @param k length of mers to be used in the background model
      @param missing_pattern_length number of missing nucleotides to be complete
      @param cur_pattern current pattern id
      @param cur_prob probability of current pattern
      @param final_probabilities array with final background probabilities of full-length patterns

  */
  void calculate_bg_probabilities(GapMask* mask, size_t* factors, float* full_mask_patterns);

  void calculate_bg_probability(float* background_model, const int alphabet_size, const int pattern_length,
                                           const int k,
                                           int missing_pattern_length, size_t pattern,
                                           float cur_prob, size_t* factors, float* final_probabilities);


  /**
      Calculation of the log(pvalues) of base patterns

      @param ltot number of possible patterns in the input
  */
  void calculate_log_pvalues(int ltot);

  /**
      Calculation of the zscore of base patterns

      @param ltot number of possible patterns in the input
  */
  void calculate_zscores(int ltot);

  /**
      filter base patterns
      - order base patterns by zscore
      - remove base patterns with edit distance one
      - remove base patterns with zscore < zscorethreshold

      @param alphabet_size the size of the alphabet
      @param zscore_threshold consider only base patterns with zscore above this threshold
      @param selected_patterns set with id's of selected base patterns
  */
  void filter_base_patterns(const int pattern_length, const int alphabet_size,
                                  const size_t number_patterns, const float zscore_threshold,
                                  const size_t count_threshold,
                                  float* pattern_zscore, std::vector<size_t>& selected_patterns);


  /**
      optimize/degenerate base patterns to IUPAC patterns

      @param selected_base_patterns set with id's of base patterns to be processed
      @param best_iupac_patterns result vector with pointers to objects of iupac patterns
  */
  void optimize_iupac_patterns(std::vector<size_t>& selected_base_patterns,
                               std::vector<IUPACPattern*>& best_iupac_patterns);

  /**
      filter iupac patterns
      - remove iupac patterns starting with N
      - remove iupac patterns with less than 4 more informative positions than non-informative positions (Ns)

      @param iupac_patterns with pointers; pointers to filtered iupac's will be removed
  */
  void filter_iupac_patterns(std::vector<IUPACPattern*>& iupac_patterns);

  /**
      merge overlapping iupac patterns

      @param iupac_patterns with pointers; merged IUPACs replace their descendants
  */
  void merge_iupac_patterns(const size_t pattern_length, const float bit_factor_merge_threshold,
                            BackgroundModel* bg, std::vector<IUPACPattern*>& iupac_patterns);

  /**
      optimize pwms in iupac patterns with respect to the number of pattern occurrences

      @param iupac_patterns with pointers; merged IUPACs replace their descendants
      @param saturation_factor a scaling factor for the em
      @param min_em_threshold the em stops if the summed up difference between the previous and current pwm are below this threshold
      @param max_iterations maximal number of iterations with the em per pattern
  */
  void em_optimize_pwms(std::vector<IUPACPattern*>& iupac_patterns,
                                  const float saturation_factor,
                                  const float min_em_threshold, const int max_iterations);

  /**
      calculates in an recursive fashion odds of probabilities needed for the em

      @param curr_pattern recursively built pattern till the pattern is complete
      @param curr_prob current probability of the pattern
      @param curr_length current length of the pattern
      @param pwm the pwm of a pattern that is optimized by the em
      @param prob_odds float array contains in the end for each base pattern the log probs
  */
  void init_prob_odds(const size_t pattern_length,
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
