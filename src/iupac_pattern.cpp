/*
 * iupac_pattern.cpp
 *
 *  Created on: Jan 30, 2017
 *      Author: mmeier
 */

#include <cstdio>
#include <iostream>
#include <cmath>
#include <limits>
#include "shared/Alphabet.h"
#include "iupac_pattern.h"
#include "iupac_alphabet.h"
#include "base_pattern.h"
#include "helper-inl.h"
#include "utils.h"

size_t* IUPACPattern::iupac_factor = nullptr;
float* IUPACPattern::log_bonferroni = nullptr;
float** IUPACPattern::iupac_profile = nullptr;

constexpr float MIXIN_FACTOR = 0.2;
constexpr float MIXIN_BIAS = 0.7;

IUPACPattern::IUPACPattern(size_t iupac_pattern, size_t pattern_length){
  this->pattern = iupac_pattern;
  this->pattern_length = pattern_length;

  this->local_n_sites = new size_t[pattern_length];

  this->optimization_bg_model_order = 0;
  this->pwm = nullptr;
  this->comp_pwm = nullptr;
  this->n_sites = 0;
  this->log_pvalue = 0;
  this->zscore = 0;
  this->bg_p = 0;
  this->expected_counts = 0;
  this->merged = false;
}

IUPACPattern::IUPACPattern(IUPACPattern* ori, float** pwm) {
  this->pattern = ori->pattern;
  this->pattern_length = ori->pattern_length;


  this->local_n_sites = new size_t[pattern_length];
  for(int i = 0; i < pattern_length; i++) {
    this->local_n_sites[i] = ori->local_n_sites[i];
  }

  this->pwm = new float*[pattern_length];
  for(int p = 0; p < pattern_length; p++) {
    this->pwm[p] = new float[4];
    for(int i = 0; i < 4; i++) {
      this->pwm[p][i] = pwm[p][i];
    }
  }
  IUPACPattern::normalize_pwm(pattern_length, this->pwm);

  this->comp_pwm = nullptr;
  calculate_comp_pwm();

  this->n_sites = ori->n_sites;
  this->log_pvalue = ori->log_pvalue;
  this->bg_p = ori->bg_p;
  this->expected_counts = ori->expected_counts;
  this->merged = ori->merged;
  this->optimization_bg_model_order = ori->optimization_bg_model_order;
}

//merge constructor
IUPACPattern::IUPACPattern(IUPACPattern* longer_pattern, IUPACPattern* shorter_pattern, bool is_comp, float* background, const int shift) {
  int offset_shorter = -1.0 * std::min(shift, 0);
  int offset_larger = std::max(shift, 0);
  int overlap = std::min(longer_pattern->get_pattern_length() - offset_larger, shorter_pattern->get_pattern_length() - offset_shorter);

  float** longer_pattern_pwm = longer_pattern->get_pwm();
  float** shorter_pattern_pwm = shorter_pattern->get_pwm();

  if(is_comp && longer_pattern->get_sites() < shorter_pattern->get_sites()) {
    longer_pattern_pwm = longer_pattern->get_comp_pwm();
  }
  else if(is_comp) {
    shorter_pattern_pwm = shorter_pattern->get_comp_pwm();
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

  this->pattern = 0;
  for(int p = 0; p < pattern_length; p++) {
    int pos_in_shorter = p - std::max(0, shift);
    int pos_in_longer = p + std::min(shift, 0);

    //only longer
    if(pos_in_shorter < 0) {
      for(int a = 0; a < 4; a++) {
        pwm[p][a] = longer_pattern_pwm[pos_in_longer][a];
      }
    }
    //only longer
    else if(pos_in_shorter >= shorter_pattern->get_pattern_length()) {
      for(int a = 0; a < 4; a++) {
        pwm[p][a] = longer_pattern_pwm[pos_in_longer][a];
      }
    }

    //only shorter
    if(pos_in_longer < 0) {
      for(int a = 0; a < 4; a++) {
        pwm[p][a] = shorter_pattern_pwm[pos_in_shorter][a];
      }
    }
    //only shorter
    else if(pos_in_longer >= longer_pattern->get_pattern_length()) {
      for(int a = 0; a < 4; a++) {
        pwm[p][a] = shorter_pattern_pwm[pos_in_shorter][a];
      }
    }

    //overlapping part
    if(pos_in_shorter >= 0 && pos_in_shorter < shorter_pattern->get_pattern_length() &&
        pos_in_longer >= 0 && pos_in_longer < longer_pattern->get_pattern_length()) {
      for(int a = 0; a < 4; a++) {
        pwm[p][a] = (shorter_pattern->local_n_sites[pos_in_shorter] * shorter_pattern_pwm[pos_in_shorter][a] +
            longer_pattern->local_n_sites[pos_in_longer] * longer_pattern_pwm[pos_in_longer][a])
            / (shorter_pattern->local_n_sites[pos_in_shorter] + longer_pattern->local_n_sites[pos_in_longer]);
      }
    }
  }

  IUPACPattern::normalize_pwm(pattern_length, pwm);
  comp_pwm = nullptr;
  this->calculate_comp_pwm();

  this->log_pvalue = calculate_merged_pvalue(longer_pattern, shorter_pattern, is_comp, background, shift);

  //TODO
  this->bg_p = 0;

  this->merged = true;
}

IUPACPattern::~IUPACPattern(){
  if(pwm) {
    for(int p = 0; p < pattern_length; p++) {
      delete[] pwm[p];
    }
    delete[] pwm;
  }

  if(comp_pwm) {
    for(int p = 0; p < pattern_length; p++) {
      delete[] comp_pwm[p];
    }
    delete[] comp_pwm;
  }

  delete[] local_n_sites;
}

void IUPACPattern::init(size_t max_pattern_length, float* bg_model) {
  //init factors to get patterns from their numerical identifiers
  iupac_factor = new size_t[max_pattern_length + 1];
  for(size_t i = 0; i < max_pattern_length + 1; i++) {
    iupac_factor[i] = pow(IUPAC_ALPHABET_SIZE, i);
  }

  log_bonferroni = new float[IUPAC_ALPHABET_SIZE];
  log_bonferroni[to_underlying(IUPAC_Alphabet::A)] = log(8);
  log_bonferroni[to_underlying(IUPAC_Alphabet::C)] = log(8);
  log_bonferroni[to_underlying(IUPAC_Alphabet::G)] = log(8);
  log_bonferroni[to_underlying(IUPAC_Alphabet::T)] = log(8);
  log_bonferroni[to_underlying(IUPAC_Alphabet::S)] = log(16);
  log_bonferroni[to_underlying(IUPAC_Alphabet::W)] = log(16);
  log_bonferroni[to_underlying(IUPAC_Alphabet::R)] = log(16);
  log_bonferroni[to_underlying(IUPAC_Alphabet::Y)] = log(16);
  log_bonferroni[to_underlying(IUPAC_Alphabet::M)] = log(24);
  log_bonferroni[to_underlying(IUPAC_Alphabet::K)] = log(24);
  log_bonferroni[to_underlying(IUPAC_Alphabet::N)] = log(6);

  initIUPACProfile(MIXIN_FACTOR, MIXIN_BIAS, bg_model);
}

void IUPACPattern::initIUPACProfile(const float mixin_factor, const float mixin_bias, float* bg_model) {
  iupac_profile = new float*[IUPAC_ALPHABET_SIZE];
  for(int i = 0; i < IUPAC_ALPHABET_SIZE; i++) {
    iupac_profile[i] = new float[4];
    for(int a = 0; a < 4; a++) {
      iupac_profile[i][a] = 0;
    }
  }


  for(int iupac_c = 0; iupac_c < IUPAC_ALPHABET_SIZE; iupac_c++) {
    std::vector<int> rep = IUPACAlphabet::get_representative_iupac_nucleotides(iupac_c);
    for(int a = 0; a < 4; a++) {
      iupac_profile[iupac_c][a] += mixin_factor * bg_model[a];

      for(auto r : rep) {
        if(a == r) {
          iupac_profile[iupac_c][a] += mixin_bias;
          break;
        }
      }
    }
  }
}

float IUPACPattern::calculate_merged_pvalue(IUPACPattern* longer_pattern, IUPACPattern* shorter_pattern, bool is_comp,
                                            float* background, const int shift) {
  float** longer_pattern_pwm = longer_pattern->get_pwm();
  float** shorter_pattern_pwm = shorter_pattern->get_pwm();

  if(is_comp && longer_pattern->get_sites() < shorter_pattern->get_sites()) {
    longer_pattern_pwm = longer_pattern->get_comp_pwm();
  }
  else {
    shorter_pattern_pwm = shorter_pattern->get_comp_pwm();
  }


  int offset_shorter = -1.0 * std::min(shift, 0);
  int offset_longer = std::max(shift, 0);
  int overlap = std::min(longer_pattern->get_pattern_length() - offset_longer, shorter_pattern->get_pattern_length() - offset_shorter);

  if(longer_pattern->log_pvalue < shorter_pattern->log_pvalue) {
    float d = 0;
    //starts with non-overlapping part
    if(offset_shorter != 0) {
      d = calculate_d_bg(shorter_pattern_pwm, background, offset_shorter, 0);
    }
    //ends with non-overlapping part
    else {
      int non_overlap_offset_start = offset_shorter + overlap;
      int non_overlap_length = shorter_pattern->get_pattern_length() - non_overlap_offset_start;
      d = calculate_d_bg(shorter_pattern_pwm, background, non_overlap_length, non_overlap_offset_start);
    }
    float d_divisor = calculate_d_bg(shorter_pattern_pwm, background, shorter_pattern->get_pattern_length());

    return longer_pattern->getLogPval() + d / d_divisor * shorter_pattern->getLogPval();
  }
  else {
    float d = 0;
    //starts with non-overlapping part
    if(offset_longer != 0) {
      d = calculate_d_bg(longer_pattern_pwm, background, offset_longer, 0);
    }
    //ends with non-overlapping part
    else {
      int non_overlap_offset_start = offset_longer + overlap;
      int non_overlap_length = longer_pattern->get_pattern_length() - non_overlap_offset_start;
      d = calculate_d_bg(longer_pattern_pwm, background, non_overlap_length, non_overlap_offset_start);
    }
    float d_divisor = calculate_d_bg(longer_pattern_pwm, background, longer_pattern->get_pattern_length());

    return shorter_pattern->getLogPval() + d / d_divisor * longer_pattern->getLogPval();
  }
}

void IUPACPattern::normalize_pwm(const int pattern_length, float** pwm) {
  const float rounding_threshold = 0.0001;
  //normalize new pwm
  for(int p = 0; p < pattern_length; p++) {
    //normalize between 0 and one -- especially necessary for em optimization
    float sum = 0;
    for(int a = 0; a < 4; a++) {
      sum += pwm[p][a];
    }
    for(int a = 0; a < 4; a++) {
      pwm[p][a] /= sum;
    }

    //get rid of too small values
    for(int a = 0; a < 4; a++) {
      if(pwm[p][a] < rounding_threshold) {
        pwm[p][a] = 0.0;
      }
    }

    //re-normalize again
    sum = 0;
    for(int a = 0; a < 4; a++) {
      sum += pwm[p][a];
    }
    for(int a = 0; a < 4; a++) {
      pwm[p][a] /= sum;
    }
  }
}


std::string IUPACPattern::toString(size_t pattern_id, size_t pattern_length) {
  std::string out = "";
  for(size_t p = 0; p < pattern_length; p++) {
    size_t residue = (pattern_id % iupac_factor[p+1]);
    out += IUPACAlphabet::getBase(residue / iupac_factor[p]);
    pattern_id -= residue;
  }
  return out;
}


int IUPACPattern::getNucleotideAtPos(const size_t pattern, const unsigned long int pos) {
  size_t residue = (pattern % IUPACPattern::iupac_factor[pos+1]);
  int c = int(residue / IUPACPattern::iupac_factor[pos]);
  return c;
}

int IUPACPattern::get_optimization_bg_model_order() {
  return optimization_bg_model_order;
}

void IUPACPattern::set_optimization_bg_model_order(int order) {
  optimization_bg_model_order = order;
}

void IUPACPattern::aggregate_attributes_from_basepatterns(BasePattern* base_pattern) {

  if(!merged) {
    base_patterns = generate_base_patterns(base_pattern, pattern);
    //works just for unmerged iupac patterns; iupac patterns need to have the same same size as the base_patterns
    int bg_order = base_pattern->getBackgroundOrder();
    float* base_background_prob = base_pattern->getBackgroundProb(bg_order);
    size_t* base_pattern_counts = base_pattern->getPatternCounter();
    float* base_pattern_exp = base_pattern->getExpectedCounts();

    // aggregate over base patterns
    float sum_backgroud_prob = 0;
    int sum_counts = 0;
    float sum_expected = 0;
    for (auto p: this->base_patterns) {
      sum_backgroud_prob += base_background_prob[p];
      sum_counts += base_pattern_counts[p];
      sum_expected += base_pattern_exp[p];
    }

    this->bg_p = sum_backgroud_prob;
    this->expected_counts = sum_expected;
    this->zscore = (sum_counts - sum_expected) / sqrt(sum_expected);
    this->n_sites = sum_counts;

    for(int p = 0; p < pattern_length; p++) {
      this->local_n_sites[p] = sum_counts;
    }

    if(n_sites == 0) {
      log_pvalue = std::numeric_limits<float>::infinity();
    } else {
      float mu = expected_counts;
      float frac = 1 - mu / (n_sites + 1);

      float log_pvalue = 0;
      if(n_sites > mu && n_sites > 5 && zscore > 2) {
        log_pvalue = n_sites * log(mu/n_sites) + n_sites - mu - 0.5 * log(6.283*n_sites*frac*frac);
      }

      for(int p = 0; p < pattern_length; p++) {
        int c = IUPACPattern::getNucleotideAtPos(pattern, p);
        log_pvalue += log_bonferroni[c];
      }

      this->log_pvalue = log_pvalue;
    }
  }

}

void IUPACPattern::calculate_pwm(BasePattern* base_pattern, const int pseudo_counts, size_t* pattern_counter, float* background_model) {
  //pwm's of merged patterns have to be initialized with the respective constructor
  if(!pwm && !merged) {
    pwm = new float*[pattern_length];
    for(int p = 0; p < pattern_length; p++) {
      pwm[p] = new float[4]; //base nucleotides ACGT
      for(int i = 0; i < 4; i++) {
        pwm[p][i] = pseudo_counts * background_model[i];
      }
    }

    for(auto base : base_patterns) {
      size_t count = pattern_counter[base];
      for(size_t p = 0; p < pattern_length; p++) {
        int c = base_pattern->getNucleotideAtPos(base, p);
        //TODO: map c to base nucleotides for MH...
        pwm[p][c] += count;
      }
    }

    for(int p = 0; p < pattern_length; p++) {
      for(int i = 0; i < 4; i++) {
        pwm[p][i] /= 1.0 * n_sites + pseudo_counts;
      }
    }

    this->calculate_comp_pwm();
  }
}

void IUPACPattern::calculate_adv_pwm(BasePattern* base_pattern, const int pseudo_counts, size_t* pattern_counter, float* background_model) {
  //pwm's of merged patterns have to be initialized with the respective constructor
  if(pwm == nullptr && !merged) {
    pwm = new float*[pattern_length];
    for(int p = 0; p < pattern_length; p++) {
      pwm[p] = new float[4]; //base nucleotides ACGT
      for(int i = 0; i < 4; i++) {
        pwm[p][i] = 0;
      }
    }

    for(size_t p = 0; p < pattern_length; p++) {
      int c = IUPACPattern::getNucleotideAtPos(pattern, p);
      size_t n_total = 0;
      size_t i_total[4];

      for(int i = 0; i < 4; i++) {
        size_t ipattern = pattern - c * IUPACPattern::iupac_factor[p] + i * IUPACPattern::iupac_factor[p];
        auto i_base_patterns = generate_base_patterns(base_pattern, ipattern);

        i_total[i] = pseudo_counts * background_model[i];
        for(auto base : i_base_patterns) {
          i_total[i] += pattern_counter[base];
        }

        n_total += i_total[i];
      }

      for(int i = 0; i < 4; i++) {
        pwm[p][i] = 1.0 * i_total[i] / n_total;
      }
    }

    this->calculate_comp_pwm();
  }
}

float IUPACPattern::calculate_d(float** p1_pwm, float** p2_pwm, const int offset1, const int offset2, const int l, const float epsilon) {
  float d = 0;
  for(int i = 0; i < l; i++) {
    for(int a = 0; a < 4; a++) {
      float mean = (p1_pwm[offset1 + i][a] + p2_pwm[offset2 + i][a] + 2 * epsilon) / 2;
      d += (p1_pwm[offset1 + i][a] + epsilon) * log2(p1_pwm[offset1 + i][a] + epsilon) + (p2_pwm[offset2 + i][a] + epsilon) * log2(p2_pwm[offset2 + i][a] + epsilon) - 2 * mean * log2(mean);
    }
  }

  return d;
}

float IUPACPattern::calculate_d_bg(float** p_pwm, float* background, const int l, const int offset, const float epsilon) {
  float d = 0;
  for(int i = 0; i < l; i++) {
    for(int a = 0; a < 4; a++) {
      float mean = (p_pwm[offset + i][a] + background[a] + 2 * epsilon) / 2;
      d += (p_pwm[offset + i][a] + epsilon) * log2(p_pwm[offset + i][a] + epsilon) + (background[a] + epsilon) * log2(background[a] + epsilon) - 2 * mean * log2(mean);
    }
  }

  return d;
}

float IUPACPattern::calculate_s(float** p1_pwm, float** p2_pwm, float* background, const int offset1, const int offset2, const int l) {
  float s = 0.5 * (calculate_d_bg(p1_pwm, background, l, offset1) +  calculate_d_bg(p2_pwm, background, l, offset2)) - calculate_d(p1_pwm, p2_pwm, offset1, offset2, l);
  return s;
}

std::tuple<float, int, bool> IUPACPattern::calculate_S(IUPACPattern* p1, IUPACPattern* p2, Strand s, float* background) {
  IUPACPattern* p_larger = p1;
  IUPACPattern* p_shorter = p2;

  if(p1->get_pattern_length() < p2->get_pattern_length()) {
    p_larger = p2;
    p_shorter = p1;
  }

  float max_s = -std::numeric_limits<float>::infinity();
  int max_shift = -255;
  bool max_comp = false;

  std::vector<bool> comp_test;
  comp_test.push_back(false);
  if(s == Strand::BOTH_STRANDS) {
    comp_test.push_back(true);
  }

  for(bool comp : comp_test) {
    for(int shift = MIN_MERGE_OVERLAP - p_shorter->get_pattern_length(); shift <= p_larger->get_pattern_length() - MIN_MERGE_OVERLAP; shift++) {
      int offset_p_shorter = -1.0 * std::min(shift, 0);
      int offset_p_larger = std::max(shift, 0);
      int overlap = std::min(p_larger->get_pattern_length() - offset_p_larger, p_shorter->get_pattern_length() - offset_p_shorter);

      float s = 0;
      if(!comp) {
        s = IUPACPattern::calculate_s(p_larger->get_pwm(), p_shorter->get_pwm(), background, offset_p_larger, offset_p_shorter, overlap);
      }
      else if(p_larger->get_sites() < p_shorter->get_sites()) {
        s = IUPACPattern::calculate_s(p_larger->get_comp_pwm(), p_shorter->get_pwm(), background, offset_p_larger, offset_p_shorter, overlap);
      }
      else { // if(p_larger->get_sites() >= p_shorter->get_sites()) {
        s = IUPACPattern::calculate_s(p_larger->get_pwm(), p_shorter->get_comp_pwm(), background, offset_p_larger, offset_p_shorter, overlap);
      }

      if(s > max_s) {
        max_s = s;
        max_shift = shift;
        max_comp = comp;
      }
    }
  }

  return std::make_tuple(max_s, max_shift, max_comp);
}


void IUPACPattern::calculate_comp_pwm() {
  if(comp_pwm == NULL) {
    comp_pwm = new float*[pattern_length];
    for(int p = 0; p < pattern_length; p++) {
      comp_pwm[p] = new float[4]; //base nucleotides ACGT
      for(int i = 0; i < 4; i++) {
        comp_pwm[p][i] = 0;
      }
    }
  }

  for(int p = 0; p < pattern_length; p++) {
    for(int i = 0; i < 4; i++) {
      comp_pwm[p][i] = pwm[pattern_length - p - 1][Alphabet::getComplementCode(i+1)-1];
    }
  }
}

void IUPACPattern::update_pwm(float** new_pwm) {
  for(int p = 0; p < pattern_length; p++) {
    for(int i = 0; i < 4; i++) {
      pwm[p][i] = new_pwm[p][i];
    }
  }

  IUPACPattern::normalize_pwm(pattern_length, pwm);
  calculate_comp_pwm();
}


float IUPACPattern::getExpCountFraction(const size_t pseudo_expected_pattern_counts) {
  return (this->expected_counts + pseudo_expected_pattern_counts) / this->n_sites;
}

float IUPACPattern::getMutualInformationScore(unsigned int n_sequences) {
  float observed_counts = this->n_sites;
  float expected_counts = this->expected_counts;

  if(observed_counts < expected_counts) {
    // these are certainly not interesting to us
    return 0;
  }

  auto MI = calculate_mutual_information_fast;
  auto H = calculate_entropy;

  float score = 0;
  for(float q: {0.5, 0.1, 0.01}) {
    score += MI(observed_counts, expected_counts, n_sequences, q) / H(q);
  }
  return -score;
}

float IUPACPattern::getLogPval() {
  return log_pvalue;
}

float IUPACPattern::getOptimizationScore(OPTIMIZATION_SCORE score_type, const size_t pseudo_expected_pattern_counts,
                                         const unsigned int n_sequences) {
  if (score_type == OPTIMIZATION_SCORE::kLogPval) {
    return getLogPval();
  } else if(score_type == OPTIMIZATION_SCORE::kExpCounts) {
    return getExpCountFraction(pseudo_expected_pattern_counts);
  }
  else if(score_type == OPTIMIZATION_SCORE::MutualInformation) {
    return getMutualInformationScore(n_sequences);
  }
  else {
    std::cerr << "Error: unknown score type!" << std::endl;
    exit(1);
  }
}

float IUPACPattern::get_bg_p() {
  return bg_p;
}

size_t IUPACPattern::get_pattern() {
  return pattern;
}

std::string IUPACPattern::get_pattern_string() {
  std::string res;

  for(int i = 0; i < pattern_length; i++) {
    float min_dist = std::numeric_limits<float>::infinity();
    int min_iupac = 0;

    for(int m = 0; m < IUPAC_ALPHABET_SIZE; m++) {
      float dist = calculate_d(pwm, iupac_profile, i, m, 1, 1E-7);
      if(dist < min_dist) {
        min_dist = dist;
        min_iupac = m;
      }
    }

    res += IUPACAlphabet::getBase(min_iupac);
  }

  return res;
}

size_t IUPACPattern::get_sites() {
  return n_sites;
}

float** IUPACPattern::get_pwm() {
  return pwm;
}

float** IUPACPattern::get_comp_pwm() {
  return comp_pwm;
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

void IUPACPattern::find_base_patterns(BasePattern* base_pattern, const size_t pattern, const size_t pattern_length, std::vector<size_t>& base_patterns) {
  std::vector<size_t> ids;
  ids.push_back(0);

  std::vector<size_t> tmp_ids;
  size_t* base_factors = base_pattern->getFactors();

  for(int p = 0; p < pattern_length; p++) {
    int c = IUPACPattern::getNucleotideAtPos(pattern, p);

    std::vector<int> representatives = IUPACAlphabet::get_representative_iupac_nucleotides(c);

    for(auto r : representatives) {
      for(auto pat : ids) {
        size_t base_pattern = pat + r * base_factors[p];
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

std::vector<size_t> IUPACPattern::generate_base_patterns(BasePattern* basepatterns, size_t iupac_pattern) {
  std::vector<size_t> kmer_stack{};
  std::vector<int> size_stack{};
  std::vector<size_t> result;

  kmer_stack.push_back(0);
  size_stack.push_back(0);
  while(kmer_stack.size() > 0) {
    auto kmer = kmer_stack.back();
    kmer_stack.pop_back();
    auto pos = size_stack.back();
    size_stack.pop_back();

    while(pos < static_cast<int>(pattern_length)) {
      int next_iupac_letter = IUPACPattern::getNucleotideAtPos(iupac_pattern, pos);
      auto letters = IUPACAlphabet::get_representative_iupac_nucleotides(next_iupac_letter);
      if (letters.size() > 1) {
        for(auto it = letters.begin() + 1; it != letters.end(); ++it) {
          kmer_stack.push_back(basepatterns->add_letter_to_the_right(kmer, pos, *it));
          size_stack.push_back(pos + 1);
        }
      }
      kmer = basepatterns->add_letter_to_the_right(kmer, pos, letters[0]);
      pos++;
    }
    result.push_back(kmer);
  }
  return result;
}

bool IUPACPattern::operator <(const IUPACPattern &rhs) const {
  return pattern < rhs.pattern;
}

float IUPACPattern::getExpectedCounts() const {
  return expected_counts;
}

float IUPACPattern::getZscore() const {
  return zscore;
}

bool sort_IUPAC_patterns(IUPACPattern* a, IUPACPattern* b) {
  return a->getLogPval() < b->getLogPval();
}

