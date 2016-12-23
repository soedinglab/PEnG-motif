#include <math.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include <set>
#include <iostream>
#include <fstream>
#include <limits>

#include "PatternCensus.h"
#include "shared/Sequence.h"
#include "iupac_alphabet.h"


int* BasePattern::factor = NULL;
int* BasePattern::rev_factor = NULL;
size_t BasePattern::pattern_length = 0;

int* IUPACPattern::iupac_factor = NULL;
float* IUPACPattern::log_bonferroni = NULL;
size_t IUPACPattern::pattern_length = 0;

void BasePattern::init(size_t pattern_length) {
  BasePattern::pattern_length = pattern_length;

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
}

std::string BasePattern::getPatternFromNumber(size_t pattern_id) {
  std::string out = "";
  for(size_t p = 0; p < pattern_length; p++) {
    int residue = (pattern_id % factor[p+1]);
    out += Alphabet::getBase((residue / factor[p]) + 1);
    pattern_id -= residue;
  }

  return out;
}

size_t BasePattern::get_rev_pattern_id(const size_t pattern_id) {
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

int BasePattern::getNucleotideAtPos(const size_t pattern, const size_t pos) {
  int residue = (pattern % factor[pos + 1]);
  int c = int(residue / factor[pos]);

  return c;
}

IUPACPattern::IUPACPattern(size_t iupac_pattern){
  pattern = iupac_pattern;
  pwm = NULL;
  n_sites = 0;
  log_pvalue = 0;
  bg_p = 0;
}

IUPACPattern::~IUPACPattern(){
  if(pwm) {
    for(int p = 0; p < pattern_length; p++) {
      delete[] pwm[p];
    }
    delete[] pwm;
  }
}

void IUPACPattern::init(size_t pattern_length) {
  IUPACPattern::pattern_length = pattern_length;

  //init factors to get patterns from their numerical identifiers
  iupac_factor = new int[pattern_length + 1];
  for(size_t i = 0; i < pattern_length + 1; i++) {
    iupac_factor[i] = pow(IUPAC_ALPHABET_SIZE, i);
  }

  float bf = log(2.0);
  log_bonferroni = new float[IUPAC_ALPHABET_SIZE];
  log_bonferroni[A] = 0.0;
  log_bonferroni[C] = 0.0;
  log_bonferroni[G] = 0.0;
  log_bonferroni[T] = 0.0;
  log_bonferroni[S] = bf;
  log_bonferroni[W] = bf;
  log_bonferroni[R] = bf;
  log_bonferroni[Y] = bf;
  log_bonferroni[M] = bf;
  log_bonferroni[K] = bf;
  log_bonferroni[N] = 0.0;
}

std::string IUPACPattern::getIUPACPatternFromNumber(size_t pattern_id) {
  std::string out = "";
  for(size_t p = 0; p < pattern_length; p++) {
    int residue = (pattern_id % iupac_factor[p+1]);
    out += IUPACAlphabet::getBase(residue / iupac_factor[p]);
    pattern_id -= residue;
  }

  return out;
}

size_t IUPACPattern::getIUPACPattern(size_t base_pattern) {
  //map pattern to basic pattern in iupac base
  //TODO: map M and H to C

  size_t iupac_pattern = 0;
  for(int p = 0; p < pattern_length; p++) {
    int c = BasePattern::getNucleotideAtPos(base_pattern, p);
    iupac_pattern += c * IUPACPattern::iupac_factor[p];
  }
  return iupac_pattern;
}

size_t IUPACPattern::getIUPACPattern(std::string base_pattern) {
  size_t iupac_pattern = 0;
  for(int p = 0; p < pattern_length; p++) {
    char nuc = base_pattern[p];
    int c = IUPACAlphabet::getCode(nuc);
    iupac_pattern += c * IUPACPattern::iupac_factor[p];
  }
  return iupac_pattern;
}


int IUPACPattern::getNucleotideAtPos(const unsigned long int pattern, const unsigned long int pos) {
  int residue = (pattern % IUPACPattern::iupac_factor[pos+1]);
  int c = int(residue / IUPACPattern::iupac_factor[pos]);
  return c;
}

void IUPACPattern::calculate_log_pvalue(const int ltot,
                                        float* base_background_prob,
                                        size_t* base_counts) {
  find_base_patterns();

  float sum_backgroud_prob = 0;
  float sum_counts = 0;
  for(auto p : base_patterns) {
    sum_backgroud_prob += base_background_prob[p];
    sum_counts += base_counts[p];
  }


  this->bg_p = sum_backgroud_prob;

  if(sum_counts == 0) {
    this->log_pvalue = std::numeric_limits<float>::infinity();
    return;
  }

  float mu = ltot * sum_backgroud_prob;
  float frac = 1 - mu / (sum_counts + 1);

  float log_pvalue = sum_counts * log(mu/sum_counts) + sum_counts - mu - 0.5 * log(6.283*sum_counts*frac*frac);

  for(int p = 0; p < pattern_length; p++) {
    int c = IUPACPattern::getNucleotideAtPos(pattern, p);
    log_pvalue += log_bonferroni[c];
  }

  this->log_pvalue = log_pvalue;
}

void IUPACPattern::count_sites(size_t* pattern_counter) {
  size_t total = 0;
  for(auto base : base_patterns) {
    total += pattern_counter[base];
  }
  this->n_sites = total;
}

void IUPACPattern::calculate_pwm(size_t* pattern_counter) {
  if(pwm == NULL) {
    pwm = new float*[pattern_length];
    for(int p = 0; p < pattern_length; p++) {
      pwm[p] = new float[4]; //base nucleotides ACGT
      for(int i = 0; i < 4; i++) {
        pwm[p][i] = 0;
      }
    }

    for(auto base : base_patterns) {
      size_t count = pattern_counter[base];
      for(size_t p = 0; p < pattern_length; p++) {
        int c = BasePattern::getNucleotideAtPos(base, p);
        //TODO: map c to base nucleotides for MH...
        pwm[p][c] += count;
      }
    }

    for(int p = 0; p < pattern_length; p++) {
      for(int i = 0; i < 4; i++) {
        pwm[p][i] /= 1.0 * n_sites;
      }
    }
  }
}

float IUPACPattern::get_log_pvalue() {
  return log_pvalue;
}

float IUPACPattern::get_bg_p() {
  return bg_p;
}


size_t IUPACPattern::get_pattern() {
  return pattern;
}

size_t IUPACPattern::get_sites() {
  return n_sites;
}

float** IUPACPattern::get_pwm() {
  return pwm;
}

std::set<size_t>& IUPACPattern::get_base_patterns() {
  return base_patterns;
}

void IUPACPattern::find_base_patterns() {
  if(base_patterns.size() > 0) {
    return;
  }

  std::set<size_t> ids;
  ids.insert(0);

  std::set<size_t> tmp_ids;

  for(int p = 0; p < pattern_length; p++) {
    int c = IUPACPattern::getNucleotideAtPos(pattern, p);

    std::vector<uint8_t> representatives = IUPACAlphabet::get_representative_iupac_nucleotides(c);

    for(auto r : representatives) {
      for(auto pat : ids) {
        size_t base_pattern = pat + r * BasePattern::factor[p];
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
  }
}

bool IUPACPattern::operator <(const IUPACPattern &rhs) const {
  return pattern < rhs.pattern;
}

PatternCensus::PatternCensus(const int pattern_length, const int k, const float zscore_threshold,
                             SequenceSet* sequence_set, BackgroundModel* bg, const char* outputFilename,
                             const char* jsonFilename) {

  BasePattern::init(pattern_length);
  IUPACPattern::init(pattern_length);
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


  count_patterns(pattern_length, Alphabet::getSize(), sequence_set);
  count_patterns_minus_strand(pattern_length, Alphabet::getSize(), pattern_counter);

  calculate_bg_probabilities(bg, alphabet_size, k);
  calculate_log_pvalues(ltot);
  calculate_zscores(ltot);

  std::set<size_t> filtered_patterns;
  filter_nearest_neighbours(alphabet_size, filtered_patterns);

  std::set<IUPACPattern*> best_iupac_patterns;
  optimize_iupac_patterns(zscore_threshold, filtered_patterns, best_iupac_patterns);

  filter_iupac_patterns(best_iupac_patterns);

  for(auto pattern : best_iupac_patterns) {
    pattern->find_base_patterns();
    pattern->count_sites(pattern_counter);
    pattern->calculate_pwm(pattern_counter);
  }

  std::string version_number("1.0.0");
  if(outputFilename != NULL) {
    std::string outfile(outputFilename);
    printShortMeme(best_iupac_patterns,
                   outfile,
                   version_number,
                   bg);
  }

  if(jsonFilename != NULL) {
    std::string json_outfile(jsonFilename);
    printJson(best_iupac_patterns,
                   json_outfile,
                   version_number,
                   bg);
  }

//  std::vector<std::string> test;
//  test.push_back(std::string("NNNNNNNN"));
//
//  for(auto pattern : test) {
//    IUPACPattern* p = new IUPACPattern(IUPACPattern::getIUPACPattern(pattern));
//    p->calculate_log_pvalue(ltot, pattern_bg_probabilities, pattern_counter);
//    p->find_base_patterns();
//    p->count_sites(pattern_counter);
//    p->calculate_pwm(pattern_counter);
//
//    std::cout << pattern << "\t" << p->get_sites() << "\t" << p->get_bg_p() << "\t" << p->get_log_pvalue() << std::endl;
//    std::set<size_t>& base_patterns = p->get_base_patterns();
//    for(auto base : base_patterns) {
//      std::cout << "\t" << BasePattern::getPatternFromNumber(base) << "\t" << pattern_counter[base] << "\t" << pattern_bg_probabilities[base] << std::endl;
//    }
//  }
}

PatternCensus::~PatternCensus() {
  delete[] pattern_counter;
  delete[] pattern_bg_probabilities;
  delete[] pattern_logp;
  delete[] pattern_zscore;
}

size_t PatternCensus::get_bg_id(const size_t pattern, const int curr_pattern_length, const int k) {
  size_t k_mer_pattern = 0;
  for(int i = curr_pattern_length - k - 1; i < curr_pattern_length; i++) {
    int c = BasePattern::getNucleotideAtPos(pattern, i);
    k_mer_pattern += c * BasePattern::factor[curr_pattern_length - i - 1];
  }
  return k_mer_pattern;
}

void PatternCensus::calculate_log_pvalues(int ltot) {
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
    float cur_prob = 1.0;
    for(int k_prime = 0; k_prime <= k; k_prime++) {
      cur_prob *= background_model[get_bg_id(pattern, k_prime+1, k_prime)];
    }
    calculate_bg_probability(background_model, alphabet_size, k, pattern_length - k - 1, pattern, cur_prob, this->pattern_bg_probabilities);
  }
}

void PatternCensus::calculate_bg_probability(float* background_model, const int alphabet_size,
                                         const int k,
                                         int missing_pattern_length, size_t pattern,
                                         float cur_prob, float* final_probabilities) {
  missing_pattern_length--;
  for(int c = 0; c < Alphabet::getSize(); c++) {
    size_t extended_pattern = pattern + c * BasePattern::factor[pattern_length - missing_pattern_length - 1];
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
          id += BasePattern::factor[p] * (seq[i+p] - 1); // -1 since the alphabet starts at 1
        }
        //could not match initial pattern without non-alphabetic character
        if(shifted) {
          i++;
          continue;
        }
      }


      //add new character to pattern
      // seq[i + pattern_length - 1] - 1; -1 since the alphabet starts at 1
      id += (seq[i + pattern_length - 1] - 1) * BasePattern::factor[pattern_length - 1];
      assert(id < number_patterns);

      size_t rev_id = BasePattern::get_rev_pattern_id(id);

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
    size_t rev_pattern_id = BasePattern::get_rev_pattern_id(i);
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
    size_t rev_pattern = BasePattern::get_rev_pattern_id(pattern);
    if(not seen_array[pattern] and not seen_array[rev_pattern]) {
      selected_patterns.insert(pattern);
      seen_array[pattern] = true;

      //iterate over neighbours and set to seen
      for(size_t p = 0; p < pattern_length; p++) {
        //mask nucleotide at position p of pattern
        int c = BasePattern::getNucleotideAtPos(pattern, p);
        size_t masked_neighbour = pattern - c * BasePattern::factor[p];

        //iterate over all possible nucleotides at position p
        for(int c = 0; c < alphabet_size; c++) {
          size_t neighbour = masked_neighbour + c * BasePattern::factor[p];
          seen_array[neighbour] = true;
        }
      }
    }
  }

  delete[] seen_array;
  delete[] sorted_array;
}

void PatternCensus::optimize_iupac_patterns(const float zscore_threshold,
                                            std::set<size_t>& selected_base_patterns,
                                            std::set<IUPACPattern*>& best_iupac_patterns) {
  std::set<size_t> seen;
  std::set<size_t> best;

  for(auto pattern : selected_base_patterns) {
    if(this->pattern_zscore[pattern] < zscore_threshold) {
      continue;
    }
//    std::cerr << "base_pattern: " << BasePattern::getPatternFromNumber(pattern) << std::endl;

    size_t iupac_pattern = IUPACPattern::getIUPACPattern(pattern);

    bool found_better_mutant = true;
    float best_log_pvalue = this->pattern_logp[pattern];
    IUPACPattern* best_mutant = new IUPACPattern(iupac_pattern);
    best_mutant->calculate_log_pvalue(ltot,
                                      this->pattern_bg_probabilities,
                                      this->pattern_counter);
    best_mutant->count_sites(pattern_counter);

    while(found_better_mutant) {
//      std::cerr << "\tbest_mutant: " << IUPACPattern::getIUPACPatternFromNumber(best_mutant->get_pattern()) << "\t" << best_mutant->get_sites() << "\t" << best_mutant->get_bg_p() << "\t" << best_log_pvalue << std::endl;
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
          IUPACPattern* mutated_pattern = new IUPACPattern(mutated_id);

          //calculate log_pvalue of mutated pattern
          mutated_pattern->calculate_log_pvalue(ltot,
                                                this->pattern_bg_probabilities,
                                                this->pattern_counter);
          mutated_pattern->count_sites(pattern_counter);

//          std::cerr << "\t\tmutated_pattern " << IUPACPattern::getIUPACPatternFromNumber(mutated_id) << "\t" << mutated_pattern->get_sites() << "\t" << mutated_pattern->get_bg_p() << "\t" << mutated_pattern->get_log_pvalue() << std::endl;

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
//        std::cerr << "\t" << IUPACPattern::getIUPACPatternFromNumber(best_mutant->get_pattern()) << " already seen!" << std::endl;
        found_better_mutant = false;
      }
      current_seen.erase(best_mutant->get_pattern());
      seen.insert(current_seen.begin(), current_seen.end());
    }

    if(best.count(best_mutant->get_pattern()) == 0 && seen.count(best_mutant->get_pattern()) == 0) {
//      std::cerr << "\tadded " << IUPACPattern::getIUPACPatternFromNumber(best_mutant->get_pattern()) << std::endl;
      best_iupac_patterns.insert(best_mutant);
      best.insert(best_mutant->get_pattern());
      seen.insert(best_mutant->get_pattern());
    }
    else{
//      std::cerr << "\tno new optimal pattern!" << std::endl;
      delete best_mutant;
      best_mutant = NULL;
    }
  }
}

void PatternCensus::filter_iupac_patterns(std::set<IUPACPattern*>& iupac_patterns) {
  std::vector<IUPACPattern*> sorted_iupac_patterns;
  sorted_iupac_patterns.insert(sorted_iupac_patterns.end(), iupac_patterns.begin(), iupac_patterns.end());
  std::sort(sorted_iupac_patterns.begin(), sorted_iupac_patterns.end(), sort_IUPAC_patterns);

  std::set<IUPACPattern*> deselected_patterns;

  for(auto pat : iupac_patterns) {
    size_t pattern = pat->get_pattern();
    int residue = (pattern % IUPACPattern::iupac_factor[0+1]);
    int c = int(residue / IUPACPattern::iupac_factor[0]);
    if(c == N) {
      deselected_patterns.insert(pat);
      continue;
    }

    int non_informative_positions = 0;
    for(int p = 1; p < pattern_length; p++) {
      int c = IUPACPattern::getNucleotideAtPos(pattern, p);
      if(c == N) {
        non_informative_positions += 1;
      }
    }

    if(pattern_length - non_informative_positions <= 3) {
      deselected_patterns.insert(pat);
    }
  }

  for(auto pat : deselected_patterns) {
    auto pat_it = iupac_patterns.find(pat);
    iupac_patterns.erase(pat_it);
    delete pat;
  }
}


void PatternCensus::printShortMeme(std::set<IUPACPattern*>& best_iupac_patterns,
                                   const std::string output_filename,
                                   const std::string version_number,
                                   BackgroundModel* bg_model) {

  std::vector<IUPACPattern*> sorted_iupac_patterns;
  sorted_iupac_patterns.insert(sorted_iupac_patterns.end(), best_iupac_patterns.begin(), best_iupac_patterns.end());
  std::sort(sorted_iupac_patterns.begin(), sorted_iupac_patterns.end(), sort_IUPAC_patterns);

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
      myfile << "MOTIF " << IUPACPattern::getIUPACPatternFromNumber(pattern->get_pattern()) << std::endl;
      myfile << "letter-probability matrix:" <<
          " alength= " << 4 <<
          " w= " << pattern_length <<
          " nsites= " << pattern->get_sites() <<
          " bg_prob= " << pattern->get_bg_p() <<
          " log(Pval)= " << pattern->get_log_pvalue()<< std::endl;

      float** pwm = pattern->get_pwm();
      for(size_t w = 0; w < pattern_length; w++) {
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

void PatternCensus::printJson(std::set<IUPACPattern*>& best_iupac_patterns,
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
    myfile << "\t\"pattern_length\" : " << pattern_length << "," << std::endl;

    myfile << "\t\"patterns\" : [" << std::endl;
    for(auto pattern : sorted_iupac_patterns) {
      myfile << "\t\t{" << std::endl;
      myfile << "\t\t\t\"iupac_motif\" : \"" << IUPACPattern::getIUPACPatternFromNumber(pattern->get_pattern()) << "\"," << std::endl;
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

bool sort_IUPAC_patterns(IUPACPattern* a, IUPACPattern* b) {
  return a->get_log_pvalue() < b->get_log_pvalue();
}
