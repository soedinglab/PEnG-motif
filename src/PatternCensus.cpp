#include <math.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include <set>
#include <iostream>
#include <fstream>
#include <limits>
#include <climits>

#include "PatternCensus.h"
#include "shared/Sequence.h"
#include "iupac_alphabet.h"

#ifdef OPENMP
  #include <omp.h>
#endif

size_t* BasePattern::factor = NULL;
size_t* BasePattern::rev_factor = NULL;
size_t BasePattern::pattern_length = 0;

size_t* IUPACPattern::iupac_factor = NULL;
float* IUPACPattern::log_bonferroni = NULL;

void BasePattern::init(size_t pattern_length) {
  BasePattern::pattern_length = pattern_length;

  //init factors to get patterns from their numerical identifiers
  factor = new size_t[pattern_length + 1];
  for(size_t i = 0; i < pattern_length + 1; i++) {
    factor[i] = pow(Alphabet::getSize(), i);
  }

  //init rev factors to get unique numerical identifiers for reverse patterns on minus strand
  //depends on alphabet_size and pattern_size
  rev_factor = new size_t[pattern_length];
  for(size_t i = 0; i < pattern_length; i++) {
    rev_factor[pattern_length - i - 1] = pow(Alphabet::getSize(), i);
  }
}

std::string BasePattern::getPatternFromNumber(size_t pattern_id) {
  std::string out = "";
  for(size_t p = 0; p < pattern_length; p++) {
    size_t residue = (pattern_id % factor[p+1]);
    out += Alphabet::getBase((residue / factor[p]) + 1);
    pattern_id -= residue;
  }

  return out;
}

size_t BasePattern::get_rev_pattern_id(const size_t pattern_id) {
  size_t tmp_id = pattern_id;
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
  size_t residue = (pattern % factor[pos + 1]);
  int c = int(residue / factor[pos]);

  return c;
}

//merge constructor
IUPACPattern::IUPACPattern(IUPACPattern* longer_pattern, IUPACPattern* shorter_pattern, float* background, const int shift) {
  int offset_shorter = -1.0 * std::min(shift, 0);
  int offset_larger = std::max(shift, 0);
  int overlap = std::min(longer_pattern->get_pattern_length() - offset_larger, shorter_pattern->get_pattern_length() - offset_shorter);
  std::cerr << "overlap: " << overlap << std::endl;

  //TODO not a real IUPAC pattern
  if(longer_pattern->get_log_pvalue() < shorter_pattern->get_log_pvalue()) {
    this->pattern = longer_pattern->pattern;
  }
  else {
    this->pattern = shorter_pattern->pattern;
  }
  this->pattern_length = longer_pattern->get_pattern_length() + shorter_pattern->get_pattern_length() - overlap;

  //caluclate local n sites
  this->local_n_sites = new size_t[pattern_length];
  for(int p = 0; p < pattern_length; p++) {
    this->local_n_sites[p] = 0;
  }
  for(int p = 0; p < shorter_pattern->get_pattern_length(); p++) {
    this->local_n_sites[std::max(shift, 0) + p] += shorter_pattern->local_n_sites[p];
  }
  for(int p = 0; p < longer_pattern->get_pattern_length(); p++) {
    this->local_n_sites[-std::min(shift, 0) + p] += longer_pattern->local_n_sites[p];
  }

  //calculate n sites
  int sum_sites = 0;
  for(int p = 0; p < pattern_length; p++) {
    sum_sites += this->local_n_sites[p];
  }
  this->n_sites = int(sum_sites / pattern_length);

  //init pwm
  this->pwm = new float*[this->pattern_length];
  for(int p = 0; p < pattern_length; p++) {
    this->pwm[p] = new float[4];
    for(int a = 0; a < 4; a++) {
      this->pwm[p][a] = 0;
    }
  }

  for(int p = 0; p < pattern_length; p++) {
    int pos_in_shorter = p - std::max(0, shift);
    int pos_in_longer = p + std::min(shift, 0);

    //only longer
    if(pos_in_shorter < 0) {
      for(int a = 0; a < 4; a++) {
        pwm[p][a] = longer_pattern->pwm[pos_in_longer][a];
      }
    }
    //only longer
    else if(pos_in_shorter >= shorter_pattern->get_pattern_length()) {
      for(int a = 0; a < 4; a++) {
        pwm[p][a] = longer_pattern->pwm[pos_in_longer][a];
      }
    }

    //only shorter
    if(pos_in_longer < 0) {
      for(int a = 0; a < 4; a++) {
        pwm[p][a] = shorter_pattern->pwm[pos_in_shorter][a];
      }
    }
    //only shorter
    else if(pos_in_longer >= longer_pattern->get_pattern_length()) {
      for(int a = 0; a < 4; a++) {
        pwm[p][a] = shorter_pattern->pwm[pos_in_shorter][a];
      }
    }

    //only shorter
    if(pos_in_shorter >= 0 && pos_in_shorter < shorter_pattern->get_pattern_length() &&
        pos_in_longer >= 0 && pos_in_longer < longer_pattern->get_pattern_length()) {
      for(int a = 0; a < 4; a++) {
        pwm[p][a] = (shorter_pattern->local_n_sites[pos_in_shorter] * shorter_pattern->pwm[pos_in_shorter][a] + longer_pattern->local_n_sites[pos_in_longer] * longer_pattern->pwm[pos_in_longer][a])
            / (shorter_pattern->local_n_sites[pos_in_shorter] + longer_pattern->local_n_sites[pos_in_longer]);
      }
    }
  }

  this->log_pvalue = calculate_merged_pvalue(longer_pattern, shorter_pattern, background, shift);

  //TODO
  this->bg_p = 0;
}

IUPACPattern::IUPACPattern(size_t iupac_pattern, size_t pattern_length){
  this->pattern = iupac_pattern;
  this->pattern_length = pattern_length;

  this->local_n_sites = new size_t[pattern_length];

  this->pwm = NULL;
  this->n_sites = 0;
  this->log_pvalue = 0;
  this->bg_p = 0;
}

IUPACPattern::~IUPACPattern(){
  if(pwm) {
    for(int p = 0; p < pattern_length; p++) {
      delete[] pwm[p];
    }
    delete[] pwm;
  }

  delete[] local_n_sites;
}

void IUPACPattern::init(size_t max_pattern_length) {
  //init factors to get patterns from their numerical identifiers
  iupac_factor = new size_t[max_pattern_length + 1];
  for(size_t i = 0; i < max_pattern_length + 1; i++) {
    iupac_factor[i] = pow(IUPAC_ALPHABET_SIZE, i);
  }

  float bf = log(2.0);
  log_bonferroni = new float[IUPAC_ALPHABET_SIZE];
  log_bonferroni[A] = log(8);
  log_bonferroni[C] = log(8);
  log_bonferroni[G] = log(8);
  log_bonferroni[T] = log(8);
  log_bonferroni[S] = log(16);
  log_bonferroni[W] = log(16);
  log_bonferroni[R] = log(16);
  log_bonferroni[Y] = log(16);
  log_bonferroni[M] = log(24);
  log_bonferroni[K] = log(24);
  log_bonferroni[N] = log(6);
}

float IUPACPattern::calculate_merged_pvalue(IUPACPattern* longer_pattern, IUPACPattern* shorter_pattern,
                                            float* background, const int shift) {
  int offset_shorter = -1.0 * std::min(shift, 0);
  int offset_longer = std::max(shift, 0);
  int overlap = std::min(longer_pattern->get_pattern_length() - offset_longer, shorter_pattern->get_pattern_length() - offset_shorter);

  if(longer_pattern->log_pvalue < shorter_pattern->log_pvalue) {
    float d = 0;
    //starts with non-overlapping part
    if(offset_shorter != 0) {
      d = calculate_d_bg(*shorter_pattern, background, offset_shorter, 0);
    }
    //ends with non-overlapping part
    else {
      int non_overlap_offset_start = offset_shorter + overlap;
      int non_overlap_length = shorter_pattern->get_pattern_length() - non_overlap_offset_start;
      d = calculate_d_bg(*shorter_pattern, background, non_overlap_length, non_overlap_offset_start);
    }
    float d_divisor = calculate_d_bg(*shorter_pattern, background, shorter_pattern->get_pattern_length());

    return longer_pattern->get_log_pvalue() + d / d_divisor * shorter_pattern->get_log_pvalue();
  }
  else {
    float d = 0;
    //starts with non-overlapping part
    if(offset_longer != 0) {
      d = calculate_d_bg(*longer_pattern, background, offset_longer, 0);
    }
    //ends with non-overlapping part
    else {
      int non_overlap_offset_start = offset_longer + overlap;
      int non_overlap_length = longer_pattern->get_pattern_length() - non_overlap_offset_start;
      d = calculate_d_bg(*longer_pattern, background, non_overlap_length, non_overlap_offset_start);
    }
    float d_divisor = calculate_d_bg(*longer_pattern, background, longer_pattern->get_pattern_length());

    return shorter_pattern->get_log_pvalue() + d / d_divisor * longer_pattern->get_log_pvalue();
  }
}

std::string IUPACPattern::getIUPACPatternFromNumber(size_t pattern_id, size_t pattern_length) {
  std::string out = "";
  for(size_t p = 0; p < pattern_length; p++) {
    size_t residue = (pattern_id % iupac_factor[p+1]);
    out += IUPACAlphabet::getBase(residue / iupac_factor[p]);
    pattern_id -= residue;
  }

  return out;
}

size_t IUPACPattern::getIUPACPattern(const size_t base_pattern, const size_t pattern_length) {
  //map pattern to basic pattern in iupac base
  //TODO: map M and H to C

  size_t iupac_pattern = 0;
  for(int p = 0; p < pattern_length; p++) {
    int c = BasePattern::getNucleotideAtPos(base_pattern, p);
    iupac_pattern += c * IUPACPattern::iupac_factor[p];
  }
  return iupac_pattern;
}

size_t IUPACPattern::getIUPACPattern(std::string base_pattern, const size_t pattern_length) {
  size_t iupac_pattern = 0;
  for(int p = 0; p < pattern_length; p++) {
    char nuc = base_pattern[p];
    int c = IUPACAlphabet::getCode(nuc);
    iupac_pattern += c * IUPACPattern::iupac_factor[p];
  }
  return iupac_pattern;
}


int IUPACPattern::getNucleotideAtPos(const size_t pattern, const unsigned long int pos) {
  size_t residue = (pattern % IUPACPattern::iupac_factor[pos+1]);
  int c = int(residue / IUPACPattern::iupac_factor[pos]);
  return c;
}

void IUPACPattern::calculate_log_pvalue(const int ltot,
                                        float* base_background_prob,
                                        size_t* base_counts) {
  find_base_patterns(pattern, pattern_length, base_patterns);

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

  for(int p = 0; p < pattern_length; p++) {
    this->local_n_sites[p] = total;
  }
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

void IUPACPattern::calculate_adv_pwm(size_t* pattern_counter) {
  if(pwm == NULL) {
    pwm = new float*[pattern_length];
    for(int p = 0; p < pattern_length; p++) {
      pwm[p] = new float[4]; //base nucleotides ACGT
      for(int i = 0; i < 4; i++) {
        pwm[p][i] = 0;
      }
    }

    for(size_t p = 0; p < pattern_length; p++) {
      int c = IUPACPattern::getNucleotideAtPos(pattern, p);
      size_t npattern = pattern - c * IUPACPattern::iupac_factor[p] + N * IUPACPattern::iupac_factor[p];
      std::vector<size_t> n_base_patterns;
      find_base_patterns(npattern, pattern_length, n_base_patterns);

      size_t total = 0;
      for(auto base : n_base_patterns) {
        size_t count = pattern_counter[base];
        total += count;
      }

      for(int i = 0; i < 4; i++) {
        size_t ipattern = pattern - c * IUPACPattern::iupac_factor[p] + i * IUPACPattern::iupac_factor[p];
        std::vector<size_t> i_base_patterns;
        find_base_patterns(ipattern, pattern_length, i_base_patterns);

        size_t itotal = 0;
        for(auto base : i_base_patterns) {
          itotal += pattern_counter[base];
        }

        //TODO: check pseudo value
        pwm[p][i] = 1.0 * itotal / total + 0.0000000001;
      }
    }
  }
}

float IUPACPattern::calculate_d(IUPACPattern& p1, IUPACPattern& p2, const int offset1, const int offset2, const int l) {
  float** p1_pwm = p1.get_pwm();
  float** p2_pwm = p2.get_pwm();

  float d = 0;
  for(int i = 0; i < l; i++) {
    for(int a = 0; a < 4; a++) {
      d += (p1_pwm[offset1 + i][a] - p2_pwm[offset2 + i][a]) * (log2(p1_pwm[offset1 + i][a]) - log2(p2_pwm[offset2 + i][a]));
    }
  }

  return d;
}

float IUPACPattern::calculate_d_bg(IUPACPattern& p, float* background, const int l, const int offset) {
  float** p_pwm = p.get_pwm();

  float d = 0;
  for(int i = 0; i < l; i++) {
    for(int a = 0; a < 4; a++) {
      float tmp = (p_pwm[i + offset][a] - background[a]) * (log2(p_pwm[i + offset][a]) - log2(background[a]));
      d += tmp;
    }
  }

  return d;
}

float IUPACPattern::calculate_s(IUPACPattern& p1, IUPACPattern& p2, float* background, const int offset1, const int offset2, const int l) {
  float s = 0.5 * (calculate_d_bg(p1, background, l) +  calculate_d_bg(p2, background, l)) - calculate_d(p1, p2, offset1, offset2, l);

  return s;
}

std::tuple<float, int> IUPACPattern::calculate_S(IUPACPattern* p1, IUPACPattern* p2, float* background) {
  IUPACPattern* p_larger = p1;
  IUPACPattern* p_shorter = p2;

  if(p1->get_pattern_length() < p2->get_pattern_length()) {
    p_larger = p2;
    p_shorter = p1;
  }

  float max_s = -std::numeric_limits<float>::infinity();
  int max_shift = -255;

  for(int shift = -2; shift < p_larger->get_pattern_length() - p_shorter->get_pattern_length() + 2; shift++) {
    int offset_p_shorter = -1.0 * std::min(shift, 0);
    int offset_p_larger = std::max(shift, 0);
    int overlap = std::min(p_larger->get_pattern_length() - offset_p_larger, p_shorter->get_pattern_length() - offset_p_shorter);

    //make sure that it is not an inclusion
    if(overlap >= p_shorter->get_pattern_length()) {
      continue;
    }

    float s = IUPACPattern::calculate_s(*p_larger, *p_shorter, background, offset_p_larger, offset_p_shorter, overlap);

    if(s > max_s) {
      max_s = s;
      max_shift = shift;
    }
  }

  return std::make_tuple(max_s, max_shift);
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

int IUPACPattern::get_pattern_length() {
  return pattern_length;
}

size_t* IUPACPattern::get_local_sites() {
  return local_n_sites;
}

std::vector<size_t>& IUPACPattern::get_base_patterns() {
  return base_patterns;
}

void IUPACPattern::find_base_patterns(const size_t pattern, const size_t pattern_length, std::vector<size_t>& base_patterns) {
  std::vector<size_t> ids;
  ids.push_back(0);

  std::vector<size_t> tmp_ids;

  for(int p = 0; p < pattern_length; p++) {
    int c = IUPACPattern::getNucleotideAtPos(pattern, p);

    std::vector<uint8_t> representatives = IUPACAlphabet::get_representative_iupac_nucleotides(c);

    for(auto r : representatives) {
      for(auto pat : ids) {
        size_t base_pattern = pat + r * BasePattern::factor[p];
        tmp_ids.push_back(base_pattern);
      }
    }

    ids.clear();
    for(auto pattern : tmp_ids) {
      ids.push_back(pattern);
    }

    tmp_ids.clear();
  }

  for(auto pattern : ids) {
    base_patterns.push_back(pattern);
  }
}

bool IUPACPattern::operator <(const IUPACPattern &rhs) const {
  return pattern < rhs.pattern;
}

PatternCensus::PatternCensus(const int pattern_length, const int k, const float zscore_threshold,
                             SequenceSet* sequence_set, BackgroundModel* bg, const char* outputFilename,
                             const char* jsonFilename) {

  int max_base_pattern_length = std::log(SIZE_MAX) / std::log(Alphabet::getSize()) - 1;
  int max_iupac_pattern_length = std::log(SIZE_MAX) / std::log(IUPAC_ALPHABET_SIZE) - 1;

  if(pattern_length > std::log(SIZE_MAX) / std::log(IUPAC_ALPHABET_SIZE) - 1 ||
      pattern_length > std::log(SIZE_MAX) / std::log(Alphabet::getSize()) - 1) {
    std::cerr << "Warning: pattern length too long!" << std::endl;
    std::cerr << "max pattern length: " << std::max(std::log(SIZE_MAX) / std::log(IUPAC_ALPHABET_SIZE) - 1, std::log(SIZE_MAX) / std::log(Alphabet::getSize()) - 1) << std::endl;
    exit(1);
  }

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

  std::vector<IUPACPattern*> best_iupac_patterns;
  optimize_iupac_patterns(zscore_threshold, filtered_patterns, best_iupac_patterns);

  filter_iupac_patterns(best_iupac_patterns);

  #pragma omp parallel for
  for(int i = 0; i < best_iupac_patterns.size(); i++) {
    IUPACPattern* pattern = best_iupac_patterns[i];
    pattern->count_sites(pattern_counter);
    pattern->calculate_adv_pwm(pattern_counter);
  }

  bool found_better = true;
  while(found_better) {
    found_better = false;

    float best_score = -std::numeric_limits<float>::infinity();
    int best_i = 0;
    int best_j = 0;
    int best_shift = 0;
    for(int i = 0; i < best_iupac_patterns.size(); i++) {
      for(int j = i+1; j < best_iupac_patterns.size(); j++) {
        IUPACPattern* p1 = best_iupac_patterns[i];
        IUPACPattern* p2 = best_iupac_patterns[j];
        auto res = IUPACPattern::calculate_S(p1, p2, bg->getV()[0]);

        if(std::get<0>(res) > best_score) {
          best_i = i;
          best_j = j;
          best_score = std::get<0>(res);
          best_shift = std::get<1>(res);
        }
      }
    }

    if(best_score > pattern_length * 0.75) {
      IUPACPattern* merged_pattern = NULL;
      if(best_iupac_patterns[best_i]->get_pattern_length() < best_iupac_patterns[best_j]->get_pattern_length()) {
        merged_pattern = new IUPACPattern(best_iupac_patterns[best_j], best_iupac_patterns[best_i],
                                                        bg->getV()[0], best_shift);
      }
      else {
        merged_pattern = new IUPACPattern(best_iupac_patterns[best_i], best_iupac_patterns[best_j],
                                                        bg->getV()[0], best_shift);
      }

      std::cerr << "Merging patterns: " << std::endl;
      std::cerr << "i: " << best_i << "\tj: " << best_j << std::endl;
      std::cerr << "shift: " << best_shift << std::endl;
      std::cerr << "\tMerging score S: " << best_score << std::endl;

      int curr_pattern_length = best_iupac_patterns[best_i]->get_pattern_length();
      float** curr_pwm = best_iupac_patterns[best_i]->get_pwm();
      size_t* curr_counts = best_iupac_patterns[best_i]->get_local_sites();
      float curr_pval = best_iupac_patterns[best_i]->get_log_pvalue();

      std::cerr << "\tpattern 1 length: " << curr_pattern_length << std::endl;
      std::cerr << "\tpattern 1 log(pval): " << curr_pval << std::endl;

      std::cerr << "\tpattern 1 local counts + pwm:" << std::endl;
      for(int p = 0; p < curr_pattern_length; p++) {
        std::cerr << "\t\t" << curr_counts[p] << "\t";
        for(int a = 0; a < 4; a++) {
          std::cerr << curr_pwm[p][a] << "\t";
        }
        std::cerr << std::endl;
      }

      curr_pattern_length = best_iupac_patterns[best_j]->get_pattern_length();
      curr_pwm = best_iupac_patterns[best_j]->get_pwm();
      curr_counts = best_iupac_patterns[best_j]->get_local_sites();
      curr_pval = best_iupac_patterns[best_j]->get_log_pvalue();

      std::cerr << "\tpattern 2 length: " << curr_pattern_length << std::endl;
      std::cerr << "\tpattern 2 log(pval): " << curr_pval << std::endl;

      std::cerr << "\tpattern 2 local counts + pwm:" << std::endl;
      for(int p = 0; p < curr_pattern_length; p++) {
        std::cerr << "\t\t" << curr_counts[p] << "\t";
        for(int a = 0; a < 4; a++) {
          std::cerr << curr_pwm[p][a] << "\t";
        }
        std::cerr << std::endl;
      }

      curr_pattern_length = merged_pattern->get_pattern_length();
      curr_pwm = merged_pattern->get_pwm();
      curr_counts = merged_pattern->get_local_sites();
      curr_pval = merged_pattern->get_log_pvalue();

      std::cerr << "\tmerged pattern length: " << curr_pattern_length << std::endl;
      std::cerr << "\tmerged pattern log(pval): " << curr_pval << std::endl;

      std::cerr << "\tmerged pattern local counts + pwm:" << std::endl;
      for(int p = 0; p < curr_pattern_length; p++) {
        std::cerr << "\t\t" << curr_counts[p] << "\t";
        for(int a = 0; a < 4; a++) {
          std::cerr << curr_pwm[p][a] << "\t";
        }
        std::cerr << std::endl;
      }

      delete best_iupac_patterns[best_j];
      delete best_iupac_patterns[best_i];

      best_iupac_patterns.erase(best_iupac_patterns.begin() + best_j);
      best_iupac_patterns.erase(best_iupac_patterns.begin() + best_i);

      best_iupac_patterns.push_back(merged_pattern);

      found_better = true;
    }
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

void PatternCensus::calculate_zscores(int ltot) {
  #pragma omp parallel for schedule(static)
  for(size_t pattern = 0; pattern < this->number_patterns; pattern++) {
    this->pattern_zscore[pattern] = this->pattern_counter[pattern] - ltot * pattern_bg_probabilities[pattern];
    this->pattern_zscore[pattern] *= this->pattern_zscore[pattern] / (ltot * pattern_bg_probabilities[pattern]);
  }
}

void PatternCensus::calculate_bg_probabilities(BackgroundModel* model, const int alphabet_size, const int k) {
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

  std::map<size_t, int> pattern_plus_positions;
  std::map<size_t, int> pattern_minus_positions;

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
                                            std::vector<IUPACPattern*>& best_iupac_patterns) {
  std::set<size_t> seen;
  std::set<size_t> best;

  for(auto pattern : selected_base_patterns) {
    if(this->pattern_zscore[pattern] < zscore_threshold) {
      continue;
    }
//    std::cerr << "base_pattern: " << BasePattern::getPatternFromNumber(pattern) << std::endl;

    size_t iupac_pattern = IUPACPattern::getIUPACPattern(pattern, pattern_length);

    bool found_better_mutant = true;
    float best_log_pvalue = this->pattern_logp[pattern];
    IUPACPattern* best_mutant = new IUPACPattern(iupac_pattern, pattern_length);
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
          IUPACPattern* mutated_pattern = new IUPACPattern(mutated_id, pattern_length);

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
      best_iupac_patterns.push_back(best_mutant);
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

void PatternCensus::filter_iupac_patterns(std::vector<IUPACPattern*>& iupac_patterns) {
  std::vector<IUPACPattern*> sorted_iupac_patterns;
  sorted_iupac_patterns.insert(sorted_iupac_patterns.end(), iupac_patterns.begin(), iupac_patterns.end());
  std::sort(sorted_iupac_patterns.begin(), sorted_iupac_patterns.end(), sort_IUPAC_patterns);

  std::set<IUPACPattern*> deselected_patterns;

  #pragma omp parallel for
  for(int i = 0; i < iupac_patterns.size(); i++) {
    IUPACPattern* pat = iupac_patterns[i];
    size_t pattern = pat->get_pattern();

    int c = IUPACPattern::getNucleotideAtPos(pattern, 0);
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
    auto pat_it = std::find(iupac_patterns.begin(), iupac_patterns.end(), pat);
    iupac_patterns.erase(pat_it);
    delete pat;
  }
}


void PatternCensus::printShortMeme(std::vector<IUPACPattern*>& best_iupac_patterns,
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
      myfile << "MOTIF " << IUPACPattern::getIUPACPatternFromNumber(pattern->get_pattern(), pattern_length) << std::endl;
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

void PatternCensus::printJson(std::vector<IUPACPattern*>& best_iupac_patterns,
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
      myfile << "\t\t\t\"iupac_motif\" : \"" << IUPACPattern::getIUPACPatternFromNumber(pattern->get_pattern(), pattern_length) << "\"," << std::endl;
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

bool sort_IUPAC_patterns(IUPACPattern* a, IUPACPattern* b) {
  return a->get_log_pvalue() < b->get_log_pvalue();
}
