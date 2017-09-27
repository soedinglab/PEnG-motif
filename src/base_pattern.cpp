/*
 * base_pattern.cpp
 *
 *  Created on: Jan 30, 2017
 *      Author: mmeier
 */

#include "shared/Alphabet.h"
#include <cmath>
#include <assert.h>
#include <map>
#include "iupac_pattern.h"
#include "utils.h"
#include <memory>
#include "base_pattern.h"

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

  n_sequences = sequence_set->getN();

  this->k = k;
  this->max_k = std::max(k, max_k);
  this->pattern_bg_probabilities = new float*[max_k+1];
  for(int background = 0; background <= max_k; background++) {
    this->pattern_bg_probabilities[background] = new float[number_patterns];
    calculate_bg_probabilities(bg, alphabet_size, background, pattern_bg_probabilities[background]);
  }

  expected_counts = new float[number_patterns];
  ltot = 0;
  if(this->strand == Strand::BOTH_STRANDS) {
    count_patterns(sequence_set);
    calculate_expected_counts();
  } else {
    count_patterns_single_strand(sequence_set);
    calculate_expected_counts_single_stranded();
  }

  calculate_log_pvalues();
  calculate_zscores();
}

BasePattern::~BasePattern() {
  delete[] pattern_counter;

  for(int i = 0; i <= max_k; i++) {
    delete[] pattern_bg_probabilities[i];
  }
  delete[] pattern_bg_probabilities;

  delete[] pattern_logp;
  delete[] pattern_zscore;
  delete[] expected_counts;
  //do not delete bg_model
}

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

float* BasePattern::getBackgroundProb() {
  return pattern_bg_probabilities[k];
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

float BasePattern::getExpCountFraction(const size_t pattern, const size_t pseudo_expected_pattern_counts) {
  return (this->expected_counts[pattern] + pseudo_expected_pattern_counts) / pattern_counter[pattern];
}

float BasePattern::getMutualInformationScore(const size_t pattern) {
  float expected_counts = this->expected_counts[pattern];
  unsigned int observed_counts = this->pattern_counter[pattern];

  auto MI = calculate_mutual_information_fast;
  auto H = calculate_entropy;

  float score = 0;
  for(float q: {0.5, 0.1, 0.01}) {
    score += MI(observed_counts, expected_counts, n_sequences, q) / H(q);
  }
  return -score;
}

float BasePattern::getLogPval(size_t pattern) {
  return pattern_logp[pattern];
}

size_t BasePattern::getLtot() {
  return ltot;
}

float BasePattern::getOptimizationScore(const OPTIMIZATION_SCORE score_type, const size_t pattern, const size_t pseudo_expected_pattern_counts) {
  if(score_type == OPTIMIZATION_SCORE::kLogPval) {
    return getLogPval(pattern);
  }
  else if(score_type == OPTIMIZATION_SCORE::kExpCounts) {
    return getExpCountFraction(pattern, pseudo_expected_pattern_counts);
  }
  else if(score_type == OPTIMIZATION_SCORE::MutualInformation) {
    return getMutualInformationScore(pattern);
  }
  else {
    std::cerr << "Error: unknown score type!" << std::endl;
    exit(1);
  }
}

int BasePattern::getBackgroundOrder() const{
  return k;
}


void BasePattern::calculate_log_pvalues() {
  #pragma omp parallel for schedule(static)
  for(size_t pattern = 0; pattern < this->number_patterns; pattern++) {
    if(this->pattern_counter[pattern] == 0) {
      this->pattern_logp[pattern] = std::numeric_limits<float>::infinity();
    }
    else {
      float mu = expected_counts[pattern];
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

void BasePattern::calculate_zscores() {
  #pragma omp parallel for schedule(static)
  for(size_t pattern = 0; pattern < this->number_patterns; pattern++) {
    this->pattern_zscore[pattern] = (this->pattern_counter[pattern] -
        expected_counts[pattern]) / sqrt(expected_counts[pattern]);
  }
}

void BasePattern::calculate_expected_counts() {
  #pragma omp parallel for schedule(static)
  for(size_t pattern = 0; pattern < this->number_patterns; pattern++) {
    if(pattern == BasePattern::getRevCompId(pattern)) {
      expected_counts[pattern] = pattern_bg_probabilities[k][pattern] * ltot;
    } else {
      expected_counts[pattern] = pattern_bg_probabilities[k][pattern] * 2 * ltot;
    }
  }
}

void BasePattern::calculate_expected_counts_single_stranded() {
  #pragma omp parallel for schedule(static)
  for(size_t pattern = 0; pattern < this->number_patterns; pattern++) {
      expected_counts[pattern] = pattern_bg_probabilities[k][pattern] * ltot;
  }
}

void BasePattern::calculate_bg_probabilities(BackgroundModel* model, const int alphabet_size, const int k, float* pattern_bg_probs) {
  float* background_model = model->getV()[k];
  size_t nr_initial_mers = pow(alphabet_size, k+1);

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t pattern = 0; pattern < nr_initial_mers; pattern++) {
    float joint_prob = 1.0;
    for(int k_prime = 0; k_prime <= k; k_prime++) {
      joint_prob *= model->getV()[k_prime][get_bg_id(pattern, k_prime+1, k_prime)];
    }
    calculate_bg_probability(background_model, alphabet_size, k, pattern_length - k - 1, pattern, joint_prob, pattern_bg_probs);
  }
}

void BasePattern::calculate_bg_probability(float* background_model, const int alphabet_size,
                                         const int k,
                                         int remaining_shifts, size_t pattern,
                                         float cur_prob, float* final_probabilities) {
  remaining_shifts--;
  size_t* base_factors = BasePattern::getFactors();
  for(int c = 0; c < Alphabet::getSize(); c++) {
    size_t extended_pattern = pattern + c * base_factors[pattern_length - remaining_shifts - 1];
    size_t kmer_id = get_bg_id(extended_pattern, pattern_length - remaining_shifts, k);

    float extended_pattern_prob = cur_prob * background_model[kmer_id];
    if(remaining_shifts == 0) {
      final_probabilities[extended_pattern] = extended_pattern_prob;
    }
    else{
      calculate_bg_probability(background_model, alphabet_size, k, remaining_shifts,
                           extended_pattern, extended_pattern_prob, final_probabilities);
    }
  }
}


//count all possible patterns of a certain length within the alphabet on the given sequences (just one strand)
//shift over characters not in the alphabet (seq[i] == 0)

void BasePattern::count_patterns(SequenceSet* sequence_set) {
  std::vector<Sequence*> sequences = sequence_set->getSequences();
  size_t* base_factors = BasePattern::getFactors();
  unsigned int* last_match_pos = new unsigned int[base_factors[pattern_length]]{};
  unsigned int j = pattern_length; // current position of W-mer in concatenated input sequences (padded with W positions in front)

  // Loop over sequences in input set
  for(size_t s = 0; s < sequences.size(); s++) {
    uint8_t* seq = sequences[s]->getSequence();
    std::unique_ptr<uint8_t[]> seq_rev = sequences[s]->createReverseComplement();
    int length = sequences[s]->getL();
    size_t id, idrev; // numerical id of the current pattern and of its reverse complement
    int p, i, irev;

    // Loop over kmer start positions in current sequence
    for (i = 0; i < length; i++, j++) {

      // Recompute pattern id from zeroth to one before last position
      for(p = 0, id = 0; p < pattern_length && seq[i] > 0 && i < length; p++, i++, j++) {
        id += base_factors[p] * (seq[i] - 1); // -1 since the alphabet starts at 1
      }
      if (p < pattern_length) continue; // for loop was terminated because seq[i] == 0 or i == length => go to next i
      idrev = BasePattern::getRevCompId(id);
      // At this position, i points to character directly right of current pattern id => L - 1 - i in rev complement seq
      irev = length - 1 - i; // next character **left** of current pattern idrev in rev complement sequence

      // Loop over pattern end positions in current sequence...
      for (; ; i++, j++, irev--) {

        // Has last match of pattern **or its reverse complement** occurred at least pattern_length positions ago?
        size_t id_min = std::min(id, idrev );
        if (last_match_pos[id_min] + pattern_length <= j) {
          pattern_counter[id_min]++; // one occurrence counted
          last_match_pos[id_min] = j; // last match position recorded
        }
        ltot++;

        if (i >= length || seq[i] == 0) break;  // exit loop if end of seq is reached or an N character is encountered

        // Compute next id and idrev
        id /= alphabet_size;
        id += (seq[i] - 1) * base_factors[pattern_length - 1];
        idrev %= base_factors[pattern_length - 1];
        idrev *= alphabet_size;
        idrev += seq_rev[irev] - 1;

        assert(id < number_patterns);
      }
      i++; j++; // if we terminated loop due to seq[i] == 0 we need to advance to next character
    }
    j += pattern_length; // match in a new sequence cannot overlap with match in previous sequence
  }
  delete[] last_match_pos;
}

void BasePattern::count_patterns_single_strand(SequenceSet* sequence_set) {
  std::vector<Sequence*> sequences = sequence_set->getSequences();
  size_t* base_factors = BasePattern::getFactors();
  unsigned int* last_match_pos = new unsigned int[base_factors[pattern_length]]{};
  unsigned int j = pattern_length; // current position of W-mer in concatenated input sequences (padded with W positions in front)

  // Loop over sequences in input set
  for(size_t s = 0; s < sequences.size(); s++) {
    uint8_t* seq = sequences[s]->getSequence();
    int length = sequences[s]->getL();
    size_t id; // numerical id of the current pattern and of its reverse complement
    int p, i;

    // Loop over kmer start positions in current sequence
    for (i = 0; i < length; i++, j++) {

      // Recompute pattern id from zeroth to one before last position
      for(p = 0, id = 0; p < pattern_length && seq[i] > 0 && i < length; p++, i++, j++) {
        id += base_factors[p] * (seq[i] - 1); // -1 since the alphabet starts at 1
      }
      if ( p < pattern_length) continue; // for loop was terminated because seq[i] == 0 or i == length => go to next i
      // At this position, i points to character directly right of current pattern id

      // Loop over pattern end positions in current sequence...
      for (; ; i++, j++) {

        // Has last match of pattern **or its reverse complement** occurred at least pattern_length positions ago?
        if (last_match_pos[id] + pattern_length <= j) {
          pattern_counter[id]++; // one occurrence counted
          last_match_pos[id] = j; // last match position recorded
        }
        ltot++;

        if (i >= length || seq[i] == 0) break;  // exit loop if end of seq is reached or an N character is encountered

        // Compute next id
        id /= alphabet_size;
        id += (seq[i] - 1) * base_factors[pattern_length - 1];
        assert(id < number_patterns);
      }
      i++; j++; // if we terminated loop due to seq[i] == 0 we need to advance to next character
    }
    j += pattern_length; // match in a new sequence cannot overlap with match in previous sequence
  }
  delete[] last_match_pos;
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

  std::cout << "basepattern" << "\t\t" << "observed" << "\t\t" << "expected" << "\t\t" << "zscore" << std::endl;
  for(auto pattern : selected_patterns) {
    std::cout << toString(pattern) << "\t\t" << pattern_counter[pattern] << "\t\t" << expected_counts[pattern]
              << "\t\t" << pattern_zscore[pattern] << std::endl;
  }

  delete[] seen_array;
  delete[] sorted_array;
}

float* BasePattern::getExpectedCounts() const {
  return expected_counts;
}




