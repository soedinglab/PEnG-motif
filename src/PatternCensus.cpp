#include <math.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include <set>

#include "PatternCensus.h"
#include "shared/Sequence.h"


PatternCensus::PatternCensus(const int pattern_length, const int k, SequenceSet* sequence_set, BackgroundModel* bg) {
  this->pattern_length = pattern_length;
  this->alphabet_size = Alphabet::getSize();

  int ltot = 0;
  std::vector<Sequence*> seqs = sequence_set->getSequences();
  for(size_t i = 0; i < seqs.size(); i++) {
    ltot += seqs[i]->getL();
  }

  //init counter for patterns
  this->number_patterns = pow(Alphabet::getSize(), pattern_length);
  this->pattern_counter = new size_t[number_patterns];
  for(size_t i = 0; i < number_patterns; i++) {
    pattern_counter[i] = 0;
  }

  this->pattern_probabilities = new float[number_patterns];
  this->pattern_logp = new float[number_patterns];
  this->pattern_zscore = new float[number_patterns];

  //init factors to get patterns from their numerical identifiers
  factor = new int[pattern_length + 1];
  for(size_t i = 0; i < pattern_length + 1; i++) {
    factor[i] = pow(Alphabet::getSize(), i);
  }

  count_patterns(pattern_length, Alphabet::getSize(), sequence_set);
  count_patterns_reverse_strand(pattern_length, Alphabet::getSize(), pattern_counter);

  calculate_bg_probabilities(bg, alphabet_size, k);
  calculate_log_pvalues(ltot);
  calculate_zscores(ltot);

  filter_nearest_neighbours(this->alphabet_size);

  std::cerr << "pattern\tcounts\tbg_prob\tlogp\tzscore" << std::endl;
  for(size_t i = 0; i < number_patterns; i++) {
    std::cerr << getPatternFromNumber(i) << "\t" << pattern_counter[i] << "\t" << pattern_probabilities[i] << "\t" << this->pattern_logp[i] << "\t" << this->pattern_zscore[i]<< std::endl;
  }
}

PatternCensus::~PatternCensus() {
  delete[] pattern_counter;
  delete[] factor;
}

std::string PatternCensus::getPatternFromNumber(size_t pattern_id) {
  std::string out = "";
  for(size_t p = 0; p < pattern_length; p++) {
    int residue = (pattern_id % factor[p+1]);
    out += Alphabet::getBase((residue / factor[p]) + 1);
    pattern_id -= residue;
  }

  return out;
}

size_t PatternCensus::get_bg_id(const size_t pattern, const int curr_pattern_length, const int k) {
  size_t k_mer_pattern = 0;
  for(int i = curr_pattern_length - k - 1; i < curr_pattern_length; i++) {
    int residue = (pattern % this->factor[i+1]);
    int c = int(residue / this->factor[i]);

    k_mer_pattern = (c+1) * this->factor[curr_pattern_length - i - 1];
  }
  return k_mer_pattern;
}

void PatternCensus::calculate_log_pvalues(int ltot) {
  for(size_t pattern = 0; pattern < this->number_patterns; pattern++) {
    float mu = ltot * this->pattern_probabilities[pattern];
    float frac = 1.0 - mu / (this->pattern_counter[pattern] + 1);
    this->pattern_logp[pattern] = this->pattern_counter[pattern] * log(mu/this->pattern_counter[pattern])
        + this->pattern_counter[pattern] - mu - 0.5 * log(6.283 * this->pattern_counter[pattern] * frac * frac);
  }
}

void PatternCensus::calculate_zscores(int ltot) {
  for(size_t pattern = 0; pattern < this->number_patterns; pattern++) {
    this->pattern_zscore[pattern] = this->pattern_counter[pattern] - ltot * pattern_probabilities[pattern];
    this->pattern_zscore[pattern] *= this->pattern_zscore[pattern] / (ltot * pattern_probabilities[pattern]);
  }
}

void PatternCensus::calculate_bg_probabilities(BackgroundModel* model, const int alphabet_size, const int k) {
  float* background_model = model->getV()[k];

  size_t nr_initial_mers = pow(alphabet_size, k+1);

  for(size_t pattern = 0; pattern < nr_initial_mers; pattern++) {
    float cur_prob = background_model[get_bg_id(pattern, k+1, k)];
    calculate_bg_probability(background_model, alphabet_size, k, pattern_length - k - 1, pattern, cur_prob, this->pattern_probabilities);
  }
}

void PatternCensus::calculate_bg_probability(float* background_model, const int alphabet_size,
                                         const int k,
                                         int missing_pattern_length, size_t pattern,
                                         float cur_prob, float* final_probabilities) {
  missing_pattern_length--;
  for(int c = 0; c < Alphabet::getSize(); c++) {
    size_t extended_pattern = pattern + c * this->factor[pattern_length - missing_pattern_length - 1];
    size_t kmer_id = get_bg_id(extended_pattern, pattern_length - missing_pattern_length, k);

    float extended_pattern_prob = cur_prob * background_model[kmer_id];
    if(missing_pattern_length == 0) {
      final_probabilities[extended_pattern] = extended_pattern_prob;
    }
    else{
      calculate_bg_probability(background_model, alphabet_size, k, missing_pattern_length,
                           extended_pattern, extended_pattern_prob, final_probabilities);
    }
  }
}


//count all possible patterns of a certain length within the alphabet on the given sequences (just one strand)
//shift over characters not in the alphabet (seq[i] == 0)
void PatternCensus::count_patterns(const int pattern_length, const int alphabet_size, SequenceSet* sequence_set) {
  std::vector<Sequence*> sequences = sequence_set->getSequences();

  for(size_t s = 0; s < sequences.size(); s++) {
    uint8_t* seq = sequences[s]->getSequence();
    int length = sequences[s]->getL();

    //sequence is too small for matches with the pattern
    if(length < pattern_length) {
      continue;
    }

    //index in sequence; start point for current pattern
    int i = 0;

    //the numerical id of the current pattern
    size_t id;

    //start point of pattern was shifted; always at the beginning and when a non-alphabetic character occurs in the sequence
    bool shifted = true;

    while(i <= length - pattern_length) {
      if(seq[i + pattern_length - 1] == 0) {
        shifted = true;
        i += pattern_length;
        continue;
      }
      if(shifted) {
        shifted = false;
        id = 0;
        //init first pattern id without character [pattern_length - 1] for sequence
        for(size_t p = 0; p < pattern_length - 1; p++) {
          //could not match initial pattern without non-alphabetic character
          if(seq[i+p] == 0) {
            shifted = true;
            break;
          }
          id += factor[p] * (seq[i+p] - 1); // -1 since the alphabet starts at 1
        }
        //could not match initial pattern without non-alphabetic character
        if(shifted) {
          i++;
          continue;
        }
      }

      //add new character to pattern
      // seq[i + pattern_length - 1] - 1; -1 since the alphabet starts at 1
      id += (seq[i + pattern_length - 1] - 1) * factor[pattern_length - 1];

      //raise counter for pattern
      assert(id < number_patterns);
      pattern_counter[id] += 1;

      //remove first character from pattern
      id -= (seq[i] - 1); // same as: -= (seq[i] - 1) * factor[0]; -1 since the alphabet starts at 1
      id /= alphabet_size;
      i++;
    }
  }
}

//add patterns from -strand; derived from the patterns in the +strand
void PatternCensus::count_patterns_reverse_strand(const int pattern_length, const int alphabet_size, size_t* pattern_counter) {
  size_t* rev_pattern_counter = new size_t[number_patterns];
  for(size_t i = 0; i < number_patterns; i++) {
    rev_pattern_counter[i] = 0;
  }

  //init rev factors to get unique numerical identifiers for reverse patterns
  //depends on alphabet_size and pattern_size
  int* rev_factor = new int[pattern_length];
  for(size_t i = 0; i < pattern_length; i++) {
    rev_factor[pattern_length - i - 1] = pow(Alphabet::getSize(), i);
  }

  uint8_t* pattern = new uint8_t[pattern_length];
  for(size_t i = 0; i < number_patterns; i++) {
    int pattern_id = i;
    //get pattern to pattern_id
    for(int p = pattern_length - 1; p >= 0; p--) {
      pattern[p] = int(pattern_id / factor[p]);
      pattern_id -= pattern[p] * factor[p];
    }

    size_t rev_pattern_id = 0;
    //get rev pattern_id to pattern
    for(size_t p = 0; p < pattern_length; p++) {
      rev_pattern_id += (Alphabet::getComplementCode(pattern[p] + 1) - 1) * rev_factor[p];
    }

    assert(rev_pattern_id < number_patterns);
    rev_pattern_counter[rev_pattern_id] = pattern_counter[i];
  }

  for(size_t i = 0; i < number_patterns; i++) {
    pattern_counter[i] += rev_pattern_counter[i];
  }

  delete[] pattern;
  delete[] rev_factor;
}

void PatternCensus::filter_nearest_neighbours(const int alphabet_size) {
  size_t* sorted_array = new size_t[this->number_patterns];
  for(size_t i = 0; i < this->number_patterns; i++) {
    sorted_array[i] = i;
  }

  bool* seen_array = new bool[this->number_patterns];
  for(size_t i = 0; i < this->number_patterns; i++) {
    seen_array[i] = false;
  }

  std::sort(this->pattern_zscore, this->pattern_zscore + this->number_patterns, sort_indices(sorted_array));
  std::set<size_t> selected;

  for(size_t i = 0; i < this->number_patterns; i++) {
    size_t pattern = sorted_array[i];
    if(not seen_array[pattern]) {
      selected.insert(pattern);
      seen_array[pattern] = true;

      //iterate over neighbours and set to seen
      for(size_t p = 0; p < pattern_length; p++) {
        //mask nucleotide at position p of pattern
        int residue = (pattern % factor[p+1]);
        size_t masked_neighbour = pattern - int(residue / factor[p]) * factor[p];

        //iterate over all possible nucleotides at position p
        for(int c = 0; c < alphabet_size; c++) {
          size_t neighbour = masked_neighbour + c * factor[p];
          seen_array[neighbour] = true;
        }
      }
    }
  }

//  std::cerr << "pattern\tcounts\tbg_prob\tlogp\tzscore" << std::endl;
//  for(auto i : selected) {
//    std::cerr << getPatternFromNumber(i) << "\t" << pattern_counter[i] << "\t" << pattern_probabilities[i] << "\t" << this->pattern_logp[i] << "\t" << this->pattern_zscore[i]<< std::endl;
//  }
}
