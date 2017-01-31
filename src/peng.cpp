#include "shared/Sequence.h"
#include <math.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include <set>
#include <iostream>
#include <fstream>
#include <limits>
#include <climits>
#include "iupac_alphabet.h"
#include "base_pattern.h"
#include "peng.h"

#ifdef OPENMP
  #include <omp.h>
#endif


Peng::Peng(const int pattern_length, const int k, SequenceSet* sequence_set, BackgroundModel* bg) {
  int max_base_pattern_length = std::log(SIZE_MAX) / std::log(Alphabet::getSize()) - 1;
  int max_iupac_pattern_length = std::log(SIZE_MAX) / std::log(IUPAC_ALPHABET_SIZE) - 1;

  if(pattern_length > std::log(SIZE_MAX) / std::log(IUPAC_ALPHABET_SIZE) - 1 ||
      pattern_length > std::log(SIZE_MAX) / std::log(Alphabet::getSize()) - 1) {
    std::cerr << "Warning: pattern length too long!" << std::endl;
    std::cerr << "max pattern length: " << std::max(std::log(SIZE_MAX) / std::log(IUPAC_ALPHABET_SIZE) - 1, std::log(SIZE_MAX) / std::log(Alphabet::getSize()) - 1) << std::endl;
    exit(1);
  }

  BasePattern::init(pattern_length);
  IUPACPattern::init(max_iupac_pattern_length);
  IUPACAlphabet::init(Alphabet::getAlphabet());

  this->pattern_length = pattern_length;
  this->alphabet_size = Alphabet::getSize();

  ltot = 0;

  std::vector<Sequence*> seqs = sequence_set->getSequences();
  for(size_t i = 0; i < seqs.size(); i++) {
    ltot += seqs[i]->getL() - pattern_length + 1;
  }
  //for the - strand
  ltot *= 2;

  //init counter for patterns
  this->number_patterns = pow(Alphabet::getSize(), pattern_length);
  this->pattern_counter = new size_t[number_patterns];
  for(size_t i = 0; i < number_patterns; i++) {
    pattern_counter[i] = 0;
  }

  this->pattern_bg_probabilities = new float[number_patterns];
  this->pattern_logp = new float[number_patterns];
  this->pattern_zscore = new float[number_patterns];

  count_patterns(pattern_length, Alphabet::getSize(), number_patterns, sequence_set, pattern_counter);
  count_patterns_minus_strand(pattern_length, Alphabet::getSize(), number_patterns, pattern_counter);

  calculate_bg_probabilities(bg, alphabet_size, k);
  calculate_log_pvalues(ltot);
  calculate_zscores(ltot);

  bg_model = bg;
}

Peng::~Peng() {
  delete[] pattern_counter;
  delete[] pattern_bg_probabilities;
  delete[] pattern_logp;
  delete[] pattern_zscore;
  //do not delete bg_model
}

size_t Peng::get_bg_id(const size_t pattern, const int curr_pattern_length, const int k) {
  size_t k_mer_pattern = 0;
  size_t* base_factors = BasePattern::getFactors();
  for(int i = curr_pattern_length - k - 1; i < curr_pattern_length; i++) {
    int c = BasePattern::getNucleotideAtPos(pattern, i);
    k_mer_pattern += c * base_factors[curr_pattern_length - i - 1];
  }
  return k_mer_pattern;
}

void Peng::calculate_log_pvalues(int ltot) {
  #pragma omp parallel for schedule(static)
  for(size_t pattern = 0; pattern < this->number_patterns; pattern++) {
    if(this->pattern_counter[pattern] == 0) {
      this->pattern_logp[pattern] = std::numeric_limits<float>::infinity();
    }
    else {
      float mu = ltot * this->pattern_bg_probabilities[pattern];
      float frac = 1.0 - mu / (this->pattern_counter[pattern] + 1);

      this->pattern_logp[pattern] = this->pattern_counter[pattern] * log(mu/this->pattern_counter[pattern])
          + this->pattern_counter[pattern] - mu - 0.5 * log(6.283 * this->pattern_counter[pattern] * frac * frac);
    }
  }
}

void Peng::calculate_zscores(int ltot) {
  #pragma omp parallel for schedule(static)
  for(size_t pattern = 0; pattern < this->number_patterns; pattern++) {
    this->pattern_zscore[pattern] = this->pattern_counter[pattern] - ltot * pattern_bg_probabilities[pattern];
    this->pattern_zscore[pattern] *= this->pattern_zscore[pattern] / (ltot * pattern_bg_probabilities[pattern]);
  }
}

void Peng::calculate_bg_probabilities(BackgroundModel* model, const int alphabet_size, const int k) {
  float* background_model = model->getV()[k];
  size_t nr_initial_mers = pow(alphabet_size, k+1);

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t pattern = 0; pattern < nr_initial_mers; pattern++) {
    float cur_prob = 1.0;
    for(int k_prime = 0; k_prime <= k; k_prime++) {
      cur_prob *= background_model[get_bg_id(pattern, k_prime+1, k_prime)];
    }
    calculate_bg_probability(background_model, alphabet_size, k, pattern_length - k - 1, pattern, cur_prob, this->pattern_bg_probabilities);
  }
}

void Peng::calculate_bg_probability(float* background_model, const int alphabet_size,
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
void Peng::count_patterns(const int pattern_length, const int alphabet_size,
                          const size_t number_patterns,
                          SequenceSet* sequence_set, size_t* pattern_counter) {
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

      size_t rev_id = BasePattern::getMinusId(id);

      if((pattern_plus_positions.find(id) == pattern_plus_positions.end() || i - pattern_plus_positions[id] >= pattern_length)
          && (pattern_minus_positions.find(id) == pattern_minus_positions.end() || i - pattern_minus_positions[id] >= pattern_length)
          && (pattern_plus_positions.find(rev_id) == pattern_plus_positions.end() || i - pattern_plus_positions[rev_id] >= pattern_length)
          && (pattern_minus_positions.find(rev_id) == pattern_minus_positions.end() || i - pattern_minus_positions[rev_id] >= pattern_length)) {
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
}

//add patterns from -strand; derived from the patterns in the +strand
void Peng::count_patterns_minus_strand(const int pattern_length, const int alphabet_size,
                                       const size_t number_patterns, size_t* pattern_counter) {
  size_t* rev_pattern_counter = new size_t[number_patterns];
  for(size_t i = 0; i < number_patterns; i++) {
    rev_pattern_counter[i] = 0;
  }

  for(size_t i = 0; i < number_patterns; i++) {
    size_t rev_pattern_id = BasePattern::getMinusId(i);
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


void Peng::filter_base_patterns(const int pattern_length, const int alphabet_size,
                                const size_t number_patterns, const float zscore_threshold,
                                float* pattern_zscore, std::set<size_t>& selected_patterns) {
  size_t* sorted_array = new size_t[number_patterns];
  for(size_t i = 0; i < number_patterns; i++) {
    sorted_array[i] = i;
  }

  bool* seen_array = new bool[number_patterns];
  for(size_t i = 0; i < number_patterns; i++) {
    seen_array[i] = false;
  }

  std::sort(sorted_array, sorted_array + number_patterns, sort_indices(pattern_zscore));
  size_t* base_factors = BasePattern::getFactors();

  for(size_t i = 0; i < number_patterns; i++) {
    size_t pattern = sorted_array[i];
    if(pattern_zscore[pattern] < zscore_threshold) {
      break;
    }

    size_t rev_pattern = BasePattern::getMinusId(pattern);
    if(not seen_array[pattern] and not seen_array[rev_pattern]) {
      selected_patterns.insert(pattern);
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

  delete[] seen_array;
  delete[] sorted_array;
}

void Peng::em_optimize_pwms(std::vector<IUPACPattern*>& best_iupac_patterns,
                                      const float saturation_factor,
                                      const float min_em_threshold, const int max_iterations) {

  #pragma omp parallel for
  for(int i = 0; i < best_iupac_patterns.size(); i++) {
    //allocate pwm's for optimization
    float** old_pwm = new float*[pattern_length];
    for(int p = 0; p < pattern_length; p++) {
      old_pwm[p] = new float[4];
    }

    float** new_pwm = new float*[pattern_length];
    for(int p = 0; p < pattern_length; p++) {
      new_pwm[p] = new float[4];
    }

    float* prob_odds = new float[number_patterns];

    //copy pwm to old_pwm
    float** ori_pwm = best_iupac_patterns[i]->get_pwm();
    for(int p = 0; p < pattern_length; p++) {
      for(int a = 0; a < 4; a++) {
        old_pwm[p][a] = ori_pwm[p][a];
      }
    }

    float change = pattern_length;
    int iteration_counter = 0;

    while(true) {
      if(change <= min_em_threshold || iteration_counter == max_iterations) {
        break;
      }

      iteration_counter++;
      //init new_pwm
      for(int p = 0; p < pattern_length; p++) {
        for(int a = 0; a < 4; a++) {
          new_pwm[p][a] = 0.0;
        }
      }

      //init prob_odds
      init_prob_odds(pattern_length, 0, 1.0, 0, old_pwm, pattern_bg_probabilities, prob_odds);

      //calculate new pwm
      for(size_t pattern = 0; pattern < number_patterns; pattern++) {
        for(int p = 0; p < pattern_length; p++) {
          int a = BasePattern::getNucleotideAtPos(pattern, p);
          new_pwm[p][a] += pattern_counter[pattern] * saturation_factor / (1 + saturation_factor / prob_odds[pattern]);
        }
      }

      //normalize new pwm
      for(int p = 0; p < pattern_length; p++) {
        float sum = 0;
        for(int a = 0; a < 4; a++) {
          sum += new_pwm[p][a];
        }
        for(int a = 0; a < 4; a++) {
          new_pwm[p][a] /= sum;
        }
      }

      //calculate change
      change = 0;
      for(int p = 0; p < pattern_length; p++) {
        for(int a = 0; a < 4; a++) {
          change += abs(new_pwm[p][a] - old_pwm[p][a]);
        }
      }

      //switch old new
      float** switcher = 0;
      switcher = old_pwm;
      old_pwm = new_pwm;
      new_pwm = switcher;
    }

    //copy new pwm to ori pwm
    for(int p = 0; p < pattern_length; p++) {
      for(int a = 0; a < 4; a++) {
        ori_pwm[p][a] = old_pwm[p][a];
      }
    }

    //de-allocate pwm's
    for(int p = 0; p < pattern_length; p++) {
      delete[] new_pwm[p];
    }
    delete[] new_pwm;

    for(int p = 0; p < pattern_length; p++) {
      delete[] old_pwm[p];
    }
    delete[] old_pwm;

    delete[] prob_odds;
  }
}

void Peng::init_prob_odds(const size_t pattern_length,
                          size_t curr_pattern, float curr_prob, int curr_length,
                          float** pwm, float* pattern_bg_probabilities, float* prob_odds) {
  if(curr_length < pattern_length) {
    size_t* factors = BasePattern::getFactors();
    for(int a = 0; a < 4; a++) {
      float next_prob = curr_prob * pwm[curr_length][a];
      size_t next_pattern = curr_pattern + a * factors[curr_length];
      int next_pattern_length = curr_length + 1;

      init_prob_odds(pattern_length, next_pattern, next_prob, next_pattern_length,
                     pwm, pattern_bg_probabilities, prob_odds);
    }
  }
  else {
    prob_odds[curr_pattern] = curr_prob / pattern_bg_probabilities[curr_pattern];
  }
}

void Peng::merge_iupac_patterns(const size_t pattern_length, const float merge_bit_factor_threshold,
                                BackgroundModel* bg, std::vector<IUPACPattern*>& best_iupac_patterns) {
  bool found_better = true;
  while(found_better) {
    found_better = false;

    //calculate all pairwise S scores between IUPACPatterns* in the vector best_iupac_patterns
    //remember the best score with the corresponding indices
    float best_score = -std::numeric_limits<float>::infinity();
    int best_i = 0;
    int best_j = 0;
    int best_shift = 0;
    for(int i = 0; i < best_iupac_patterns.size(); i++) {
      IUPACPattern* p1 = best_iupac_patterns[i];
      if(p1->get_log_pvalue() > -5) {
        continue;
      }
      for(int j = i+1; j < best_iupac_patterns.size(); j++) {
        IUPACPattern* p2 = best_iupac_patterns[j];
        if(p2->get_log_pvalue() > -5) {
          continue;
        }

        auto res = IUPACPattern::calculate_S(p1, p2, bg_model->getV()[0]);

        if(std::get<0>(res) > best_score) {
          best_i = i;
          best_j = j;
          best_score = std::get<0>(res);
          best_shift = std::get<1>(res);
        }
      }
    }

    //if best score is above threshold; merge and continue searching for appropriate merges
    if(best_score > pattern_length * merge_bit_factor_threshold) {
      IUPACPattern* merged_pattern = NULL;
      if(best_iupac_patterns[best_i]->get_pattern_length() < best_iupac_patterns[best_j]->get_pattern_length()) {
        merged_pattern = new IUPACPattern(best_iupac_patterns[best_j], best_iupac_patterns[best_i],
                                                        bg_model->getV()[0], best_shift);
      }
      else {
        merged_pattern = new IUPACPattern(best_iupac_patterns[best_i], best_iupac_patterns[best_j],
                                                        bg_model->getV()[0], best_shift);
      }

      int curr_pattern_length = best_iupac_patterns[best_i]->get_pattern_length();
      float** curr_pwm = best_iupac_patterns[best_i]->get_pwm();
      size_t* curr_counts = best_iupac_patterns[best_i]->get_local_sites();
      float curr_pval = best_iupac_patterns[best_i]->get_log_pvalue();

      curr_pattern_length = best_iupac_patterns[best_j]->get_pattern_length();
      curr_pwm = best_iupac_patterns[best_j]->get_pwm();
      curr_counts = best_iupac_patterns[best_j]->get_local_sites();
      curr_pval = best_iupac_patterns[best_j]->get_log_pvalue();

      curr_pattern_length = merged_pattern->get_pattern_length();
      curr_pwm = merged_pattern->get_pwm();
      curr_counts = merged_pattern->get_local_sites();
      curr_pval = merged_pattern->get_log_pvalue();

      delete best_iupac_patterns[best_j];
      delete best_iupac_patterns[best_i];

      best_iupac_patterns.erase(best_iupac_patterns.begin() + best_j);
      best_iupac_patterns.erase(best_iupac_patterns.begin() + best_i);

      best_iupac_patterns.push_back(merged_pattern);

      //continue searching for merges
      found_better = true;
    }
  }
}

void Peng::process(const float zscore_threshold,
                   const bool use_em, const float em_saturation_factor, const float min_em_threshold,
                   const int em_max_iterations, const float bit_factor_merge_threshold,
                   std::vector<IUPACPattern*>& best_iupac_patterns) {
  std::set<size_t> filtered_base_patterns;
  filter_base_patterns(pattern_length, Alphabet::getSize(), number_patterns, zscore_threshold, pattern_zscore, filtered_base_patterns);

  optimize_iupac_patterns(filtered_base_patterns, best_iupac_patterns);

  filter_iupac_patterns(best_iupac_patterns);

  #pragma omp parallel for
  for(int i = 0; i < best_iupac_patterns.size(); i++) {
    IUPACPattern* pattern = best_iupac_patterns[i];
    pattern->count_sites(pattern_counter);
    pattern->calculate_adv_pwm(pattern_counter, bg_model->getV()[0]);
  }

  if(use_em) {
    em_optimize_pwms(best_iupac_patterns, em_saturation_factor, min_em_threshold, em_max_iterations);
  }

  if(pattern_length >= MIN_MERGE_OVERLAP) {
    merge_iupac_patterns(pattern_length, bit_factor_merge_threshold, bg_model, best_iupac_patterns);
  }
  else {
    std::cerr << "Warning: Specified pattern length ("
        << pattern_length << ") is too low for merging!" << std::endl;
  }
}

void Peng::optimize_iupac_patterns(std::set<size_t>& selected_base_patterns,
                                            std::vector<IUPACPattern*>& best_iupac_patterns) {
  std::set<size_t> seen;
  std::set<size_t> best;

  for(auto pattern : selected_base_patterns) {
    size_t iupac_pattern = IUPACPattern::baseToId(pattern, pattern_length);

    bool found_better_mutant = true;
    float best_log_pvalue = this->pattern_logp[pattern];
    IUPACPattern* best_mutant = new IUPACPattern(iupac_pattern, pattern_length);
    best_mutant->calculate_log_pvalue(ltot,
                                      this->pattern_bg_probabilities,
                                      this->pattern_counter);
    best_mutant->count_sites(pattern_counter);

    while(found_better_mutant) {
      found_better_mutant = false;
      size_t mutant_mother = best_mutant->get_pattern();
      std::set<size_t> current_seen;

      for(int p = 0; p < pattern_length; p++) {
        //mask nucleotide at position i
        int c = IUPACPattern::getNucleotideAtPos(mutant_mother, p);
        size_t masked_mother = mutant_mother - c * IUPACPattern::iupac_factor[p];

        //replace position p with similar IUPAC nucleotides
        for(auto r : IUPACAlphabet::get_similar_iupac_nucleotides(c)) {
          size_t mutated_id = masked_mother + r * IUPACPattern::iupac_factor[p];
          IUPACPattern* mutated_pattern = new IUPACPattern(mutated_id, pattern_length);

          //calculate log_pvalue of mutated pattern
          mutated_pattern->calculate_log_pvalue(ltot,
                                                this->pattern_bg_probabilities,
                                                this->pattern_counter);
          mutated_pattern->count_sites(pattern_counter);

          //add pattern to currently seen
          current_seen.insert(mutated_pattern->get_pattern());

          if(mutated_pattern->get_log_pvalue() < best_log_pvalue) {
            delete best_mutant;
            found_better_mutant = true;
            best_log_pvalue = mutated_pattern->get_log_pvalue();
            best_mutant = mutated_pattern;
          }
          else {
            delete mutated_pattern;
          }
        }
      }

      if(seen.count(best_mutant->get_pattern()) == 1) {
        found_better_mutant = false;
      }
      current_seen.erase(best_mutant->get_pattern());
      seen.insert(current_seen.begin(), current_seen.end());
    }

    if(best.count(best_mutant->get_pattern()) == 0 && seen.count(best_mutant->get_pattern()) == 0) {
      best_iupac_patterns.push_back(best_mutant);
      best.insert(best_mutant->get_pattern());
      seen.insert(best_mutant->get_pattern());
    }
    else{
      delete best_mutant;
      best_mutant = NULL;
    }
  }
}

void Peng::filter_iupac_patterns(std::vector<IUPACPattern*>& iupac_patterns) {
  std::vector<IUPACPattern*> sorted_iupac_patterns;
  sorted_iupac_patterns.insert(sorted_iupac_patterns.end(), iupac_patterns.begin(), iupac_patterns.end());
  std::sort(sorted_iupac_patterns.begin(), sorted_iupac_patterns.end(), sort_IUPAC_patterns);

  std::set<IUPACPattern*> deselected_patterns;

  #pragma omp parallel for
  for(int i = 0; i < iupac_patterns.size(); i++) {
    IUPACPattern* pat = iupac_patterns[i];
    size_t pattern = pat->get_pattern();

    //pattern shall not start with non-informative position ('N')
    int c = IUPACPattern::getNucleotideAtPos(pattern, 0);
    if(c == N) {
      deselected_patterns.insert(pat);
      continue;
    }

    //count non-informative positions ('N')
    int non_informative_positions = 0;
    for(int p = 1; p < pattern_length; p++) {
      int c = IUPACPattern::getNucleotideAtPos(pattern, p);
      if(c == N) {
        non_informative_positions += 1;
      }
    }

    //limit fraction of non-informative positions ('N')
    if(pattern_length - non_informative_positions <= 3) {
      deselected_patterns.insert(pat);
    }
  }

  for(auto pat : deselected_patterns) {
    auto pat_it = std::find(iupac_patterns.begin(), iupac_patterns.end(), pat);
    iupac_patterns.erase(pat_it);
    delete pat;
  }
}


void Peng::printShortMeme(std::vector<IUPACPattern*>& best_iupac_patterns,
                                   const std::string output_filename,
                                   const std::string version_number,
                                   BackgroundModel* bg_model) {

  std::vector<IUPACPattern*> sorted_iupac_patterns;
  sorted_iupac_patterns.insert(sorted_iupac_patterns.end(), best_iupac_patterns.begin(), best_iupac_patterns.end());
  std::sort(sorted_iupac_patterns.begin(), sorted_iupac_patterns.end(), sort_IUPAC_patterns);

  auto res = IUPACPattern::calculate_S(sorted_iupac_patterns[0], sorted_iupac_patterns[1], bg_model->getV()[0]);
  float best_score = std::get<0>(res);
  int best_shift = std::get<1>(res);

  std::ofstream myfile (output_filename);
  if (myfile.is_open()) {
    myfile << "PEnG-motif version " << version_number << std::endl;
    myfile << std::endl;

    char* alphabet = Alphabet::getAlphabet();
    myfile << "ALPHABET= " << alphabet << std::endl;
    myfile << std::endl;

    myfile << "Background letter frequencies" << std::endl;
    float* freq_nuc = bg_model->getV()[0];

    for(size_t i = 0; i < strlen(alphabet); i++) {
      if(i != 0) {
        myfile << " ";
      }
      myfile << alphabet[i] << " " << freq_nuc[i];
    }
    myfile << std::endl;
    myfile << std::endl;

    for(auto pattern : sorted_iupac_patterns) {
      myfile << "MOTIF " << IUPACPattern::toString(pattern->get_pattern(), BasePattern::getPatternLength()) << std::endl;
      myfile << "letter-probability matrix:" <<
          " alength= " << 4 <<
          " w= " << pattern->get_pattern_length() <<
          " nsites= " << pattern->get_sites() <<
          " bg_prob= " << pattern->get_bg_p() <<
          " log(Pval)= " << pattern->get_log_pvalue()<< std::endl;

      float** pwm = pattern->get_pwm();
      for(size_t w = 0; w < pattern->get_pattern_length(); w++) {
        for(size_t a = 0; a < 4; a++) {
          if(a != 0) {
            myfile << " ";
          }
          myfile << pwm[w][a];
        }
        myfile << std::endl;
      }
      myfile << std::endl;
    }
    myfile.close();
  }
  else std::cerr << "Unable to open output file (" << output_filename << ")!";
}

void Peng::printJson(std::vector<IUPACPattern*>& best_iupac_patterns,
                                   const std::string output_filename,
                                   const std::string version_number,
                                   BackgroundModel* bg_model) {

  std::vector<IUPACPattern*> sorted_iupac_patterns;
  sorted_iupac_patterns.insert(sorted_iupac_patterns.end(), best_iupac_patterns.begin(), best_iupac_patterns.end());
  std::sort(sorted_iupac_patterns.begin(), sorted_iupac_patterns.end(), sort_IUPAC_patterns);

  std::ofstream myfile (output_filename);
  if (myfile.is_open()) {
    myfile << "{" << std::endl;
    char* alphabet = Alphabet::getAlphabet();
    myfile << "\t\"alphabet\" : \"" << alphabet << "\"," << std::endl;
    myfile << "\t\"bg\" : [";
    float* freq_nuc = bg_model->getV()[0];

    for(size_t i = 0; i < strlen(alphabet); i++) {
      myfile << freq_nuc[i];
      if(i != strlen(alphabet) - 1) {
        myfile << ", ";
      }
    }
    myfile << "]," << std::endl;

    myfile << "\t\"alphabet_length\" : " << 4 << "," << std::endl;

    myfile << "\t\"patterns\" : [" << std::endl;
    for(auto pattern : sorted_iupac_patterns) {
      myfile << "\t\t{" << std::endl;
      myfile << "\t\t\t\"iupac_motif\" : \"" << IUPACPattern::toString(pattern->get_pattern(), BasePattern::getPatternLength()) << "\"," << std::endl;
      myfile << "\t\t\t\"pattern_length\" : " << pattern->get_pattern_length() << "," << std::endl;
      myfile << "\t\t\t\"sites\" : " << pattern->get_sites() << "," << std::endl;
      myfile << "\t\t\t\"log(Pval)\" : " << pattern->get_log_pvalue() << "," << std::endl;
      myfile << "\t\t\t\"bg_prob\" : " << pattern->get_bg_p() << "," << std::endl;
      myfile << "\t\t\t\"pwm\" : [" << std::endl;
      float** pwm = pattern->get_pwm();
      for(size_t w = 0; w < pattern_length; w++) {
        myfile << "\t\t\t\t\t[";
        for(size_t a = 0; a < 4; a++) {
          myfile << pwm[w][a];
          if(a != 3) {
            myfile << ", ";
          }
          else{
            myfile << "]";
          }
        }
        if(w != pattern_length - 1) {
          myfile << ", ";
        }
        myfile << std::endl;
      }
      myfile << "\t\t\t\t]" << std::endl;
      myfile << "\t\t}";
      if(pattern != sorted_iupac_patterns[sorted_iupac_patterns.size() - 1]) {
        myfile << ",";
      }
      myfile << std::endl;
    }
    myfile << "\t]" << std::endl;
    myfile << "}" << std::endl;

    myfile.close();
  }
  else std::cerr << "Unable to open output file (" << output_filename << ")!";
}
