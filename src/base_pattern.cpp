/*
 * base_pattern.cpp
 *
 *  Created on: Jan 30, 2017
 *      Author: mmeier
 */

#include "shared/Alphabet.h"
#include <cstdio>
#include <cmath>
#include <assert.h>
#include <map>
#include "iupac_pattern.h"
#include "base_pattern.h"

void BasePattern::init(size_t pattern_length) {
  BasePattern::pattern_length = pattern_length;

  //init factors to get patterns from their numerical identifiers
  //"binary" encoding of patterns
  factor = new size_t[pattern_length + 1];
  for(size_t i = 0; i < pattern_length + 1; i++) {
    factor[i] = pow(Alphabet::getSize(), i);
  }
}

std::string BasePattern::toString(size_t pattern_id) {
  std::string out = "";
  for(size_t p = 0; p < pattern_length; p++) {
    int c = BasePattern::getNucleotideAtPos(pattern_id, p);
    //+1; shifted encoding compared to Alphabet (Alphabet encodes 'other' on position 0)
    out += Alphabet::getBase(c + 1);
  }

  return out;
}

size_t BasePattern::getRevCompId(const size_t pattern_id) {
  size_t tmp_id = pattern_id;
  size_t rev_pattern_id = 0;

  //calculate id of pattern on the minus strand
  for(int p = 0; p < pattern_length; p++) {
    int c = BasePattern::getNucleotideAtPos(pattern_id, p);

    //+1; shifted encoding compared to Alphabet (Alphabet encodes 'other' on position 0)
    //-1; to get back to our encoding (without 'other')
    //reverse pattern order factor[pattern_length - 1 - p]
    //reverse nucleotides Alphabet::getComplementCode
    rev_pattern_id += (Alphabet::getComplementCode(c + 1) - 1) * factor[pattern_length - 1 - p];
  }

  return rev_pattern_id;
}

int BasePattern::getNucleotideAtPos(const size_t pattern, const size_t pos) {
  //id: a0*|a|^0 + a1*|a|^1 + a2*|a|^2 + a3*|a|^3
  // a2 = floor((id % |a|^3) / |a|^2)
  size_t residue = (pattern % factor[pos + 1]);
  return int(residue / factor[pos]);
}

size_t* BasePattern::getFactors() {
  return BasePattern::factor;
}

size_t BasePattern::getPatternLength() {
  return BasePattern::pattern_length;
}

size_t BasePattern::getNumberPatterns() {
  return number_patterns;
}

size_t* BasePattern::getPatternCounter() {
  return pattern_counter;
}

float* BasePattern::getBackgroundProb(const int order) {
  return pattern_bg_probabilities[order];
}

size_t BasePattern::baseId2IUPACId(const size_t base_pattern) {
  //map pattern to basic pattern in iupac base
  size_t iupac_pattern = 0;
  for (int p = 0; p < pattern_length; p++) {
    int c = getNucleotideAtPos(base_pattern, p);
    iupac_pattern += c * IUPACPattern::iupac_factor[p];
  }
  return iupac_pattern;
}

float BasePattern::getLogPval(size_t pattern) {
  return pattern_logp[pattern];
}

size_t BasePattern::getLtot() {
  return ltot;
}

BasePattern::BasePattern(const size_t pattern_length, Strand s, const int k, const int max_k,
                         SequenceSet* sequence_set, BackgroundModel* bg) {
  this->pattern_length = pattern_length;
  alphabet_size = Alphabet::getSize();
  strand = s;

  this->init(pattern_length);

  //init counter for patterns
  number_patterns = pow(Alphabet::getSize(), pattern_length);
  pattern_counter = new size_t[number_patterns];
  for(size_t i = 0; i < number_patterns; i++) {
    pattern_counter[i] = 0;
  }

  pattern_logp = new float[number_patterns];
  pattern_zscore = new float[number_patterns];

  count_patterns(sequence_set);
  if(this->strand == Strand::BOTH_STRANDS) {
    count_patterns_minus_strand();
  }

  ltot = 0;
  for(size_t i = 0; i < number_patterns; i++) {
    ltot += pattern_counter[i];
  }

  this->k = k;
  this->max_k = std::max(k, max_k);
  this->pattern_bg_probabilities = new float*[max_k+1];
  for(int background = 0; background <= max_k; background++) {
    this->pattern_bg_probabilities[background] = new float[number_patterns];
    calculate_bg_probabilities(bg, alphabet_size, background, pattern_bg_probabilities[background]);
  }

  calculate_log_pvalues(ltot);
  calculate_zscores(ltot);
}

BasePattern::~BasePattern() {
  delete[] pattern_counter;

  for(int i = 0; i <= max_k; i++) {
    delete[] pattern_bg_probabilities[i];
  }
  delete[] pattern_bg_probabilities;

  delete[] pattern_logp;
  delete[] pattern_zscore;
  //do not delete bg_model
}

size_t BasePattern::get_bg_id(const size_t pattern, const int curr_pattern_length, const int k) {
  size_t k_mer_pattern = 0;
  size_t* base_factors = BasePattern::getFactors();
  for(int i = curr_pattern_length - k - 1; i < curr_pattern_length; i++) {
    int c = BasePattern::getNucleotideAtPos(pattern, i);
    k_mer_pattern += c * base_factors[curr_pattern_length - i - 1];
  }
  return k_mer_pattern;
}

void BasePattern::calculate_log_pvalues(int ltot) {
  #pragma omp parallel for schedule(static)
  for(size_t pattern = 0; pattern < this->number_patterns; pattern++) {
    if(this->pattern_counter[pattern] == 0) {
      this->pattern_logp[pattern] = std::numeric_limits<float>::infinity();
    }
    else {
      float mu = ltot * this->pattern_bg_probabilities[k][pattern];
      float frac = 1.0 - mu / (this->pattern_counter[pattern] + 1);

      if(this->pattern_counter[pattern] > mu && this->pattern_counter[pattern] > 5) {
        this->pattern_logp[pattern] = this->pattern_counter[pattern] * log(mu/this->pattern_counter[pattern])
            + this->pattern_counter[pattern] - mu - 0.5 * log(6.283 * this->pattern_counter[pattern] * frac * frac);
      }
      else {
        this->pattern_logp[pattern] = 0;
      }
    }
  }
}

void BasePattern::calculate_zscores(int ltot) {
  #pragma omp parallel for schedule(static)
  for(size_t pattern = 0; pattern < this->number_patterns; pattern++) {
    this->pattern_zscore[pattern] = (this->pattern_counter[pattern] -
        ltot * pattern_bg_probabilities[k][pattern]) / sqrt(ltot * pattern_bg_probabilities[k][pattern]);
  }
}

void BasePattern::calculate_bg_probabilities(BackgroundModel* model, const int alphabet_size, const int k, float* pattern_bg_probs) {
  float* background_model = model->getV()[k];
  size_t nr_initial_mers = pow(alphabet_size, k+1);

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t pattern = 0; pattern < nr_initial_mers; pattern++) {
    float cur_prob = 1.0;
    for(int k_prime = 0; k_prime <= k; k_prime++) {
      cur_prob *= model->getV()[k_prime][get_bg_id(pattern, k_prime+1, k_prime)];
    }
    calculate_bg_probability(background_model, alphabet_size, k, pattern_length - k - 1, pattern, cur_prob, pattern_bg_probs);
  }
}

void BasePattern::calculate_bg_probability(float* background_model, const int alphabet_size,
                                         const int k,
                                         int missing_pattern_length, size_t pattern,
                                         float cur_prob, float* final_probabilities) {
  missing_pattern_length--;
  size_t* base_factors = BasePattern::getFactors();
  for(int c = 0; c < Alphabet::getSize(); c++) {
    size_t extended_pattern = pattern + c * base_factors[pattern_length - missing_pattern_length - 1];
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
void BasePattern::count_patterns(SequenceSet* sequence_set) {
  std::vector<Sequence*> sequences = sequence_set->getSequences();

  std::map<size_t, int> pattern_plus_positions;
  std::map<size_t, int> pattern_minus_positions;
  size_t* base_factors = BasePattern::getFactors();

  for(size_t s = 0; s < sequences.size(); s++) {
    uint8_t* seq = sequences[s]->getSequence();
    int length = sequences[s]->getL();

    //sequence is too small for matches with the pattern
    if(length < pattern_length) {
      continue;
    }

    pattern_plus_positions.clear();
    pattern_minus_positions.clear();

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
          id += base_factors[p] * (seq[i+p] - 1); // -1 since the alphabet starts at 1
        }
        //could not match initial pattern without non-alphabetic character
        if(shifted) {
          i++;
          continue;
        }
      }

      //add new character to pattern
      // seq[i + pattern_length - 1] - 1; -1 since the alphabet starts at 1
      id += (seq[i + pattern_length - 1] - 1) * base_factors[pattern_length - 1];
      assert(id < number_patterns);

//      for(int j = 0; j < i; j++) {
//        std::cerr << " ";
//      }
//      std::cerr << BasePattern::toString(id) << std::endl;

      size_t rev_id = BasePattern::getRevCompId(id);

      if((pattern_plus_positions.find(id) == pattern_plus_positions.end() || i - pattern_plus_positions[id] >= pattern_length)
          && (pattern_minus_positions.find(id) == pattern_minus_positions.end() || i - pattern_minus_positions[id] >= pattern_length)
          && (pattern_plus_positions.find(rev_id) == pattern_plus_positions.end() || i - pattern_plus_positions[rev_id] >= pattern_length)
          && (pattern_minus_positions.find(rev_id) == pattern_minus_positions.end() || i - pattern_minus_positions[rev_id] >= pattern_length)) {

        //raise counter for pattern
        pattern_counter[id] += 1;
//        std::cerr << "counts: " << pattern_counter[id] << std::endl;

        pattern_plus_positions[id] = i;
        pattern_minus_positions[id] = i;
        pattern_plus_positions[rev_id] = i;
        pattern_minus_positions[rev_id] = i;
      }

      //remove first character from pattern
      id -= (seq[i] - 1); // same as: -= (seq[i] - 1) * factor[0]; -1 since the alphabet starts at 1
      id /= alphabet_size;
      i++;
    }
  }
}

//add patterns from -strand; derived from the patterns in the +strand
void BasePattern::count_patterns_minus_strand() {
  size_t* rev_pattern_counter = new size_t[number_patterns];
  for(size_t i = 0; i < number_patterns; i++) {
    rev_pattern_counter[i] = 0;
  }

  for(size_t i = 0; i < number_patterns; i++) {
    size_t rev_pattern_id = BasePattern::getRevCompId(i);
    assert(rev_pattern_id < number_patterns);

    //do not count palindromic patterns from minus strand
    if(rev_pattern_id != i) {
      rev_pattern_counter[rev_pattern_id] = pattern_counter[i];
    }
  }

  for(size_t i = 0; i < number_patterns; i++) {
    pattern_counter[i] += rev_pattern_counter[i];
  }

  delete[] rev_pattern_counter;
}


void BasePattern::filter_base_patterns(const float zscore_threshold,
                                       const size_t count_threshold,
                                       std::vector<size_t>& selected_patterns) {

  bool* seen_array = new bool[number_patterns];
  for(size_t i = 0; i < number_patterns; i++) {
    seen_array[i] = false;
  }

  size_t* sorted_array = new size_t[number_patterns];
  for(size_t i = 0; i < number_patterns; i++) {
    sorted_array[i] = i;
  }
  std::sort(sorted_array, sorted_array + number_patterns, sort_indices(pattern_zscore));
  size_t* base_factors = BasePattern::getFactors();

  for(size_t i = 0; i < number_patterns; i++) {
    size_t pattern = sorted_array[i];
    if(pattern_zscore[pattern] < zscore_threshold) {
      break;
    }
    if(pattern_counter[pattern] < count_threshold) {
      continue;
    }

    size_t rev_pattern = BasePattern::getRevCompId(pattern);
    if(not seen_array[pattern] and not seen_array[rev_pattern]) {
      selected_patterns.push_back(pattern);
      seen_array[pattern] = true;

      //iterate over neighbours and set to seen
      for(size_t p = 0; p < pattern_length; p++) {
        //mask nucleotide at position p of pattern
        int c = BasePattern::getNucleotideAtPos(pattern, p);
        size_t masked_neighbour = pattern - c * base_factors[p];

        //iterate over all possible nucleotides at position p
        for(int c = 0; c < alphabet_size; c++) {
          size_t neighbour = masked_neighbour + c * base_factors[p];
          seen_array[neighbour] = true;
        }
      }
    }
  }

  for(auto pattern : selected_patterns) {
    std::cout << "selected base pattern: " << toString(pattern) << "\t" << pattern_counter[pattern] << "\t" << pattern_zscore[pattern] << "\t" << pattern_logp[pattern] << std::endl;
  }

  delete[] seen_array;
  delete[] sorted_array;
}




