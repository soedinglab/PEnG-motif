#ifndef PATTERN_CENSUS_H
#define PATTERN_CENSUS_H

#include "shared/SequenceSet.h"
#include "shared/Alphabet.h"
#include "shared/BackgroundModel.h"
#include <set>
#include <array>
#include <tuple>
#include "iupac_pattern.h"

class Peng{
 public:
  Peng(const int pattern_length, const int k,
                SequenceSet* sequence_set, BackgroundModel* bg);
  ~Peng();
  void process(std::vector<IUPACPattern*>& best_iupac_patterns, const float zscore_threshold);

  void printShortMeme(std::vector<IUPACPattern*>& best_iupac_patterns,
                      const std::string output_filename,
                      const std::string version_number,
                      BackgroundModel* bg_model);

  void printJson(std::vector<IUPACPattern*>& best_iupac_patterns,
                 const std::string output_filename,
                 const std::string version_number,
                 BackgroundModel* bg_model);

 private:
  void count_patterns(const int pattern_length, const int alphabet_size, SequenceSet* sequence_set);
  void count_patterns_minus_strand(const int pattern_length, const int alphabet_size, size_t* pattern_counter);

  /**
      Helper function for the background model
      Extracts a k_mer id from a pattern to

      @param pattern is a base pattern id
      @param curr_pattern_length is the current position in the pattern
      @param k is the length of the kmer
      @return the id of the k_mer ending at curr_pattern_pos
  */
  size_t get_bg_id(const size_t pattern, const int curr_pattern_length, const int k);

  /**
      Calculates the background probabilities of base patterns

      @param model the background model
      @param alphabet_size the size of the alphabet
      @param k length of mers to be used in the background model
  */
  void calculate_bg_probabilities(BackgroundModel* model, const int alphabet_size, const int k);

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
  void calculate_bg_probability(float* background_model, const int alphabet_size,
                            const int k, int missing_pattern_length, size_t cur_pattern,
                            float cur_prob, float* final_probabilities);

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
  void filter_nearest_neighbours(const int alphabet_size, const float zscore_threshold, std::set<size_t>& selected_patterns);

  /**
      optimize/degenerate base patterns to IUPAC patterns

      @param selected_base_patterns set with id's of base patterns to be processed
      @param best_iupac_patterns result vector with pointers to objects of iupac patterns
  */
  void optimize_iupac_patterns(std::set<size_t>& selected_base_patterns,
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
  void merge_iupac_patterns(std::vector<IUPACPattern*>& iupac_patterns);

  BackgroundModel* bg_model;
  size_t* pattern_counter;
  float* pattern_bg_probabilities;
  float* pattern_logp;
  float* pattern_zscore;

  size_t number_patterns;
  int pattern_length;
  int alphabet_size;
  size_t ltot;
};

class sort_indices {
   private:
     float* mparr;
   public:
     sort_indices(float* parr) : mparr(parr) {}
     bool operator()(const size_t i, const size_t j) const { return mparr[i] > mparr[j]; }
};

#endif /* UTIL_H_ */
