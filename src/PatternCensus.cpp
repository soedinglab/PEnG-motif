#include <math.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include <set>
#include <iostream>
#include <fstream>

#include "PatternCensus.h"
#include "shared/Sequence.h"
#include "iupac_alphabet.h"


PatternCensus::PatternCensus(const int pattern_length, const int k, SequenceSet* sequence_set, BackgroundModel* bg, const char* outputFilename) {
  this->pattern_length = pattern_length;
  this->alphabet_size = Alphabet::getSize();

  ltot = 0;
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

  this->pattern_bg_probabilities = new float[number_patterns];
  this->pattern_logp = new float[number_patterns];
  this->pattern_zscore = new float[number_patterns];

  //init factors to get patterns from their numerical identifiers
  factor = new int[pattern_length + 1];
  for(size_t i = 0; i < pattern_length + 1; i++) {
    factor[i] = pow(Alphabet::getSize(), i);
  }

  //init rev factors to get unique numerical identifiers for reverse patterns on minus strand
  //depends on alphabet_size and pattern_size
  rev_factor = new int[pattern_length];
  for(size_t i = 0; i < pattern_length; i++) {
    rev_factor[pattern_length - i - 1] = pow(Alphabet::getSize(), i);
  }

  //init factors to get patterns from their numerical identifiers
  iupac_factor = new int[pattern_length + 1];
  for(size_t i = 0; i < pattern_length + 1; i++) {
    iupac_factor[i] = pow(IUPAC_ALPHABET_SIZE, i);
  }

  float bf = log(2.0);
  log_bonferroni = new float[IUPAC_ALPHABET_SIZE];
  log_bonferroni[0] = 0.0;
  log_bonferroni[1] = 0.0;
  log_bonferroni[2] = 0.0;
  log_bonferroni[3] = 0.0;
  log_bonferroni[4] = bf;
  log_bonferroni[5] = bf;
  log_bonferroni[6] = bf;
  log_bonferroni[7] = bf;
  log_bonferroni[8] = bf;
  log_bonferroni[9] = bf;
  log_bonferroni[10] = 0.0;

  IUPACAlphabet::init(Alphabet::getAlphabet());

  count_patterns(pattern_length, Alphabet::getSize(), sequence_set);
  count_patterns_minus_strand(pattern_length, Alphabet::getSize(), pattern_counter);

  calculate_bg_probabilities(bg, alphabet_size, k);
  calculate_log_pvalues(ltot);
  calculate_zscores(ltot);

  std::set<size_t> filtered_patterns;
  filter_nearest_neighbours(this->alphabet_size, filtered_patterns);

  std::set<size_t> best_iupac_patterns;
  optimize_iupac_patterns(pattern_length, filtered_patterns, best_iupac_patterns);

  filter_iupac_patterns(best_iupac_patterns, pattern_length, iupac_factor);

  std::ofstream myfile (outputFilename);
  if (myfile.is_open()) {
    for(auto pattern : best_iupac_patterns) {
      myfile << getIUPACPatternFromNumber(pattern) << std::endl;
    }
    myfile.close();
  }
  else std::cerr << "Unable to open output file (" << outputFilename << ")!";

  //get_pwm_for_iupac_patterns(best_iupac_patterns, pattern_length, pattern_counter);
}

PatternCensus::~PatternCensus() {
  delete[] pattern_counter;
  delete[] pattern_bg_probabilities;
  delete[] pattern_logp;
  delete[] pattern_zscore;

  delete[] factor;
  delete[] rev_factor;
  delete[] iupac_factor;
  delete[] log_bonferroni;
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

std::string PatternCensus::getIUPACPatternFromNumber(size_t pattern_id) {
  std::string out = "";
  for(size_t p = 0; p < pattern_length; p++) {
    int residue = (pattern_id % iupac_factor[p+1]);
    out += IUPACAlphabet::getBase(residue / iupac_factor[p]);
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
    float mu = ltot * this->pattern_bg_probabilities[pattern];
    float frac = 1.0 - mu / (this->pattern_counter[pattern] + 1);
    this->pattern_logp[pattern] = this->pattern_counter[pattern] * log(mu/this->pattern_counter[pattern])
        + this->pattern_counter[pattern] - mu - 0.5 * log(6.283 * this->pattern_counter[pattern] * frac * frac);
  }
}

void PatternCensus::calculate_zscores(int ltot) {
  for(size_t pattern = 0; pattern < this->number_patterns; pattern++) {
    this->pattern_zscore[pattern] = this->pattern_counter[pattern] - ltot * pattern_bg_probabilities[pattern];
    this->pattern_zscore[pattern] *= this->pattern_zscore[pattern] / (ltot * pattern_bg_probabilities[pattern]);
  }
}

void PatternCensus::calculate_bg_probabilities(BackgroundModel* model, const int alphabet_size, const int k) {
  float* background_model = model->getV()[k];

  size_t nr_initial_mers = pow(alphabet_size, k+1);

  for(size_t pattern = 0; pattern < nr_initial_mers; pattern++) {
    float cur_prob = background_model[get_bg_id(pattern, k+1, k)];
    calculate_bg_probability(background_model, alphabet_size, k, pattern_length - k - 1, pattern, cur_prob, this->pattern_bg_probabilities);
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

  int* pattern_plus_positions = new int[number_patterns];
  int* pattern_minus_positions = new int[number_patterns];

  for(size_t s = 0; s < sequences.size(); s++) {
    uint8_t* seq = sequences[s]->getSequence();
    int length = sequences[s]->getL();


    //sequence is too small for matches with the pattern
    if(length < pattern_length) {
      continue;
    }

    //init last pattern positions with -pattern_length
    for(int i = 0; i < number_patterns; i++) {
      pattern_plus_positions[i] = -pattern_length;
      pattern_minus_positions[i] = -pattern_length;
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
      assert(id < number_patterns);

      size_t rev_id = get_rev_pattern_id(id, pattern_length, factor, rev_factor);

      if(i - pattern_plus_positions[id] >= pattern_length
          && i - pattern_minus_positions[id] >= pattern_length
          && i - pattern_plus_positions[rev_id] >= pattern_length
          && i - pattern_minus_positions[rev_id] >= pattern_length) {
        //raise counter for pattern
        pattern_counter[id] += 1;

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

  delete[] pattern_plus_positions;
  delete[] pattern_minus_positions;
}

//add patterns from -strand; derived from the patterns in the +strand
void PatternCensus::count_patterns_minus_strand(const int pattern_length, const int alphabet_size, size_t* pattern_counter) {
  size_t* rev_pattern_counter = new size_t[number_patterns];
  for(size_t i = 0; i < number_patterns; i++) {
    rev_pattern_counter[i] = 0;
  }

  for(size_t i = 0; i < number_patterns; i++) {
    size_t rev_pattern_id = get_rev_pattern_id(i, pattern_length, factor, rev_factor);
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

size_t PatternCensus::get_rev_pattern_id(const size_t pattern_id, const int pattern_length,
                                         const int* factor, const int* rev_factor) {
  int tmp_id = pattern_id;
  size_t rev_pattern_id = 0;
  //get pattern to pattern_id
  for(int p = pattern_length - 1; p >= 0; p--) {
    int c = int(tmp_id / factor[p]);
    rev_pattern_id += (Alphabet::getComplementCode(c + 1) - 1) * rev_factor[p];
    tmp_id -= c * factor[p];
  }

  return rev_pattern_id;
}

void PatternCensus::filter_nearest_neighbours(const int alphabet_size,
                                              std::set<size_t>& selected_patterns) {
  size_t* sorted_array = new size_t[this->number_patterns];
  for(size_t i = 0; i < this->number_patterns; i++) {
    sorted_array[i] = i;
  }

  bool* seen_array = new bool[this->number_patterns];
  for(size_t i = 0; i < this->number_patterns; i++) {
    seen_array[i] = false;
  }

  std::sort(sorted_array, sorted_array + this->number_patterns, sort_indices(pattern_zscore));

  for(size_t i = 0; i < this->number_patterns; i++) {
    size_t pattern = sorted_array[i];
    size_t rev_pattern = get_rev_pattern_id(pattern, pattern_length, factor, rev_factor);
    if(not seen_array[pattern] and not seen_array[rev_pattern]) {
      selected_patterns.insert(pattern);
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
//  for(auto i : selected_patterns) {
//    std::cerr << getPatternFromNumber(i) << "\t" << pattern_counter[i] << "\t" << pattern_probabilities[i] << "\t" << this->pattern_logp[i] << "\t" << this->pattern_zscore[i]<< std::endl;
//  }

  delete[] seen_array;
  delete[] sorted_array;
}

void PatternCensus::optimize_iupac_patterns(const int pattern_length,
                                            std::set<size_t>& selected_patterns,
                                            std::set<size_t>& best_iupac_patterns) {
  std::set<size_t> seen;

  for(auto pattern : selected_patterns) {
    if(this->pattern_zscore[pattern] < 1000) {
      continue;
    }
    //map pattern to basic pattern in iupac base
    //TODO: map M and H to C
    size_t iupac_pattern = 0;
    for(int p = 0; p < pattern_length; p++) {
      int residue = (pattern % factor[p+1]);
      char c = int(residue / factor[p]);

      iupac_pattern += c * iupac_factor[p];
    }

    bool found_better_mutant = true;
    float best_log_pvalue = this->pattern_logp[pattern];
    size_t best_mutant = iupac_pattern;

    while(found_better_mutant) {
      found_better_mutant = false;
      size_t mutant_mother = best_mutant;
      std::set<size_t> current_seen;

      for(int p = 0; p < pattern_length; p++) {
        //mask nucleotide at position i
        int residue = (mutant_mother % iupac_factor[p+1]);
        char c = int(residue / iupac_factor[p]);
        size_t masked_mother = mutant_mother - c * iupac_factor[p];

        //replace position p with similar IUPAC nucleotides
        for(auto r : IUPACAlphabet::get_similar_iupac_nucleotides(c)) {
          size_t mutated_pattern = masked_mother + r * iupac_factor[p];

          //calculate zscore of mutated pattern
          float log_pvalue = calculate_logpvalue_of_iupac_pattern(mutated_pattern, ltot,
                                                           this->pattern_bg_probabilities,
                                                           this->pattern_counter);

          if(log_pvalue < best_log_pvalue) {
            found_better_mutant = true;
            best_log_pvalue = log_pvalue;
            best_mutant = mutated_pattern;
          }

          //add pattern to currently seen
          current_seen.insert(mutated_pattern);
        }
      }

      if(seen.count(best_mutant) == 1) {
        found_better_mutant = false;
      }
      seen.insert(current_seen.begin(), current_seen.end());
    }

    best_iupac_patterns.insert(best_mutant);
  }
}

float PatternCensus::calculate_logpvalue_of_iupac_pattern(size_t mutated_pattern, const int ltot,
                                                         float* base_background_prob, size_t* base_counts) {
  std::set<size_t> base_patterns;
  find_base_patterns(mutated_pattern, pattern_length, factor, iupac_factor, base_patterns);

  float sum_backgroud_prob = 0;
  float sum_counts = 0;
  for(auto p : base_patterns) {
    sum_backgroud_prob += base_background_prob[p];
    sum_counts += base_counts[p];
  }

  float mu = ltot * sum_backgroud_prob;
  float frac = 1 - mu / (sum_counts + 1);

  float log_pvalue = sum_counts * log(mu/sum_counts) + sum_counts - mu - 0.5 * log(6.283*sum_counts*frac*frac);

  for(int p = 0; p < pattern_length; p++) {
    int residue = (mutated_pattern % iupac_factor[p+1]);
    int c = int(residue / iupac_factor[p]);

    log_pvalue += log_bonferroni[c];
  }

  return log_pvalue;
}

float PatternCensus::find_base_patterns(const size_t mutated_pattern, const int pattern_length,
                                        int* factor, int* iupac_factor, std::set<size_t>& base_patterns) {
  std::set<size_t> ids;
  ids.insert(0);

  std::set<size_t> tmp_ids;

//  std::cerr << getIUPACPatternFromNumber(mutated_pattern) << std::endl;

  for(int p = 0; p < pattern_length; p++) {
    int residue = (mutated_pattern % iupac_factor[p+1]);
    int c = int(residue / iupac_factor[p]);

    std::vector<uint8_t> representatives = IUPACAlphabet::get_representative_iupac_nucleotides(c);

    for(auto r : representatives) {
      for(auto pattern : ids) {
        size_t base_pattern = pattern + r * factor[p];
        tmp_ids.insert(base_pattern);
      }
    }

    ids.clear();
    for(auto pattern : tmp_ids) {
      ids.insert(pattern);
    }

    tmp_ids.clear();
  }

  for(auto pattern : ids) {
    base_patterns.insert(pattern);
//    std::cerr << "\tbase pattern:" << getPatternFromNumber(pattern) << std::endl;
  }
}

void PatternCensus::filter_iupac_patterns(std::set<size_t>& iupac_patterns, const int pattern_length, int* iupac_factor) {
  std::set<size_t> selected_patterns;
  for(auto pattern : iupac_patterns) {
    int residue = (pattern % iupac_factor[0+1]);
    int c = int(residue / iupac_factor[0]);
    if(c == N) {
      continue;
    }

    int non_informative_positions = 0;
    for(int p = 1; p < pattern_length; p++) {
      int residue = (pattern % iupac_factor[p+1]);
      int c = int(residue / iupac_factor[p]);
      if(c == N) {
        non_informative_positions += 1;
      }
    }

    if(pattern_length - non_informative_positions <= 3) {
      continue;
    }

    selected_patterns.insert(pattern);
  }

  iupac_patterns.clear();
  for(auto pattern : selected_patterns) {
    iupac_patterns.insert(pattern);
  }
}

void PatternCensus::get_pwm_for_iupac_patterns(std::set<size_t>& best_iupac_patterns, const int pattern_length, size_t* pattern_counter) {
  for(auto pattern : best_iupac_patterns) {
    std::set<size_t> base_patterns;
    find_base_patterns(pattern, pattern_length, factor, iupac_factor, base_patterns);

    float** pwm = new float*[pattern_length];
    for(int p = 0; p < pattern_length; p++) {
      pwm[p] = new float[4]; //base nucleotides ACGT
      for(int i = 0; i < 4; i++) {
        pwm[p][i] = 0;
      }
    }

    size_t total = 0;
    for(auto base : base_patterns) {
      size_t count = pattern_counter[base];
      total += count;
      for(size_t p = 0; p < pattern_length; p++) {
        int residue = (base % factor[p+1]);
        int c = int(residue / factor[p]);

        //TODO: map c to base nucleotides for MH...
        pwm[p][c] += count;
      }
    }

    float pvalue = calculate_logpvalue_of_iupac_pattern(pattern, ltot, pattern_bg_probabilities, pattern_counter);

    std::cerr << getIUPACPatternFromNumber(pattern) << "\tlog(pvalue): " << pvalue << std::endl;
    for(int p = 0; p < pattern_length; p++) {
      std::cerr << p;
      for(int i = 0; i < 4; i++) {
        pwm[p][i] /= 1.0 * total;
        std::cerr << "\t" << pwm[p][i];
      }
      std::cerr << std::endl;
    }

    for(int p = 0; p < pattern_length; p++) {
      delete[] pwm[p];
    }
    delete[] pwm;
  }
}
