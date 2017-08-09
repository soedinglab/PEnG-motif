#include "shared/Sequence.h"
#include <vector>
#include <cassert>
#include <algorithm>
#include <set>
#include <iostream>
#include <cmath>
#include <fstream>
#include <limits>
#include <climits>
#include "iupac_alphabet.h"
#include "base_pattern.h"
#include "peng.h"
#include "helper-inl.h"

#ifdef OPENMP
  #include <omp.h>
#endif


Peng::Peng(Strand s, const int k, const int max_opt_k,
           SequenceSet* sequence_set, BackgroundModel* bg) {

  int max_base_pattern_length = log(SIZE_MAX) / log(Alphabet::getSize()) - 1;
  int max_iupac_pattern_length = log(SIZE_MAX) / log(IUPAC_ALPHABET_SIZE) - 1;

  this->alphabet_size = Alphabet::getSize();
  this->k = k;
  this->max_k = std::max(k, max_opt_k);
  this->strand = s;

  IUPACAlphabet::init(Alphabet::getAlphabet());

  //IUPACAlphabet needs to be initialized before this
  IUPACPattern::init(max_iupac_pattern_length, bg->getV()[0]);

  bg_model = bg;
  this->sequence_set = sequence_set;
}

Peng::~Peng() {
  //do not delete bg_model, sequence_set
}


void Peng::em_optimize_pwms(std::vector<IUPACPattern*>& best_iupac_patterns,
                            BasePattern* base_patterns,
                            const float saturation_factor,
                            const float min_em_threshold,
                            const int max_iterations,
                            float* bg_probabilities,
                            std::vector<IUPACPattern*>& optimized_iupac_patterns) {

  const size_t pattern_length = base_patterns->getPatternLength();
  const size_t number_patterns = base_patterns->getNumberPatterns();
  size_t* pattern_counter = base_patterns->getPatternCounter();

  //TODO threads
  int threads = 1;
  //allocate pwm's for optimization
  float*** threads_old_pwm = new float**[threads];
  for(int thread = 0; thread < threads; thread++) {
    threads_old_pwm[thread] = new float*[pattern_length];
    for(int p = 0; p < pattern_length; p++) {
      threads_old_pwm[thread][p] = new float[4];
    }
  }

  float*** threads_new_pwm = new float**[threads];
  for(int thread = 0; thread < threads; thread++) {
    threads_new_pwm[thread] = new float*[pattern_length];
    for(int p = 0; p < pattern_length; p++) {
      threads_new_pwm[thread][p] = new float[4];
    }
  }

  float** threads_prob_odds = new float*[threads];
  for(int thread = 0; thread < threads; thread++) {
    threads_prob_odds[thread] = new float[number_patterns];
  }

  for(int i = 0; i < best_iupac_patterns.size(); i++) {
    size_t ori_pattern = best_iupac_patterns[i]->get_pattern();

    //TODO threads
    int thread_number = 0;
    float** new_pwm = threads_new_pwm[thread_number];
    float** old_pwm = threads_old_pwm[thread_number];
    float* prob_odds= threads_prob_odds[thread_number];

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
      calculate_prob_odds(pattern_length, 0, 1.0, 0, old_pwm, bg_probabilities, base_patterns->getFactors(), prob_odds);

      //calculate new pwm
      for(size_t pattern = 0; pattern < number_patterns; pattern++) {
        for(int p = 0; p < pattern_length; p++) {
          int a = base_patterns->getNucleotideAtPos(pattern, p);
          new_pwm[p][a] += pattern_counter[pattern] * saturation_factor / (1 + saturation_factor / prob_odds[pattern]);//saturation_factor / (1 + saturation_factor / prob_odds[pattern]);// prob_odds[pattern]
        }
      }

      IUPACPattern::normalize_pwm(pattern_length, new_pwm);

      //calculate change
      change = 0;
      for(int p = 0; p < pattern_length; p++) {
        for(int a = 0; a < 4; a++) {
          change += std::abs(new_pwm[p][a] - old_pwm[p][a]);
        }
      }

      //switch old new
      float** switcher = 0;
      switcher = old_pwm;
      old_pwm = new_pwm;
      new_pwm = switcher;
    }

    IUPACPattern* optimized_pattern = new IUPACPattern(best_iupac_patterns[i], old_pwm);
    optimized_iupac_patterns.push_back(optimized_pattern);
    std::cout << "em: " << IUPACPattern::toString(ori_pattern, pattern_length) << " -> " << optimized_pattern->get_pattern_string() << std::endl;
  }

  //de-allocate pwm's
  for(int thread = 0; thread < threads; thread++) {
    for(int p = 0; p < pattern_length; p++) {
      delete[] threads_new_pwm[thread][p];
    }
    delete[] threads_new_pwm[thread];
  }
  delete[] threads_new_pwm;

  for(int thread = 0; thread < threads; thread++) {
    for(int p = 0; p < pattern_length; p++) {
      delete[] threads_old_pwm[thread][p];
    }
    delete[] threads_old_pwm[thread];
  }
  delete[] threads_old_pwm;

  for(int thread = 0; thread < threads; thread++) {
    delete [] threads_prob_odds[thread];

  }
  delete[] threads_prob_odds;
}

void Peng::calculate_prob_odds(const size_t pattern_length,
                          size_t curr_pattern, float curr_prob, int curr_length,
                          float** pwm, float* pattern_bg_probabilities, size_t* factors, float* prob_odds) {
  if(curr_length < pattern_length) {
    for(int a = 0; a < 4; a++) {
      float next_prob = curr_prob * pwm[curr_length][a];
      size_t next_pattern = curr_pattern + a * factors[curr_length];
      int next_pattern_length = curr_length + 1;

      calculate_prob_odds(pattern_length, next_pattern, next_prob, next_pattern_length,
                     pwm, pattern_bg_probabilities, factors, prob_odds);
    }
  }
  else {
    prob_odds[curr_pattern] = curr_prob / pattern_bg_probabilities[curr_pattern];
  }
}

void Peng::filter_redundancy(const float merge_bit_factor_threshold, std::vector<IUPACPattern*>& iupac_patterns) {
  std::sort(iupac_patterns.begin(), iupac_patterns.end(), sort_IUPAC_patterns);

  std::set<size_t, std::greater<size_t>> deselected;

  for(size_t i = 0; i < iupac_patterns.size(); i++) {
    if(deselected.find(i) == deselected.end()) {
      for(size_t j = i+1; j < iupac_patterns.size(); j++) {
        if(deselected.find(j) == deselected.end() &&
          iupac_patterns[i]->get_pattern_length() == iupac_patterns[j]->get_pattern_length()) {

          //TODO depends on strand...
          float similarity_score = IUPACPattern::calculate_s(iupac_patterns[i]->get_pwm(),
                                                           iupac_patterns[j]->get_pwm(),
                                                           bg_model->getV()[0], 0, 0,
                                                           iupac_patterns[i]->get_pattern_length());
          float comp_similarity_score = IUPACPattern::calculate_s(iupac_patterns[i]->get_comp_pwm(),
                                                                         iupac_patterns[j]->get_pwm(),
                                                                         bg_model->getV()[0], 0, 0,
                                                                         iupac_patterns[i]->get_pattern_length());
          float threshold = merge_bit_factor_threshold * iupac_patterns[i]->get_pattern_length();
          if(similarity_score > threshold || comp_similarity_score > threshold) {
            deselected.insert(j);
            break;
          }
        }
      }
    }
  }

  for (std::set<size_t>::iterator it=deselected.begin(); it != deselected.end(); ++it) {
    size_t index = *it;
    delete iupac_patterns[index];
    iupac_patterns.erase(iupac_patterns.begin() + index);
  }
}

void Peng::merge_iupac_patterns(const size_t pattern_length, const float merge_bit_factor_threshold,
                                BackgroundModel* bg, std::vector<IUPACPattern*>& iupac_patterns) {
  bool found_better = true;
  //TODO not necessary to recalculate the whole grid each time...
  while(found_better) {
    found_better = false;

    //calculate all pairwise S scores between IUPACPatterns* in the vector best_iupac_patterns
    //remember the best score with the corresponding indices
    float best_score = -std::numeric_limits<float>::infinity();
    int best_i = 0;
    int best_j = 0;
    int best_shift = 0;
    bool best_comp = false;
    for(int i = 0; i < iupac_patterns.size(); i++) {
      IUPACPattern* p1 = iupac_patterns[i];
      if(p1->getLogPval() > -5) {
        continue;
      }
      for(int j = i+1; j < iupac_patterns.size(); j++) {
        IUPACPattern* p2 = iupac_patterns[j];
        if(p2->getLogPval() > -5) {
          continue;
        }

        auto res = IUPACPattern::calculate_S(p1, p2, strand, bg_model->getV()[0]);

        if(std::get<0>(res) > best_score) {
          best_i = i;
          best_j = j;
          best_score = std::get<0>(res);
          best_shift = std::get<1>(res);
          best_comp = std::get<2>(res);
        }
      }
    }

    //if best score is above threshold; merge and continue searching for appropriate merges
    if(best_score > pattern_length * merge_bit_factor_threshold) {
      IUPACPattern* merged_pattern = nullptr;
      if(iupac_patterns[best_i]->get_pattern_length() < iupac_patterns[best_j]->get_pattern_length()) {
        merged_pattern = new IUPACPattern(iupac_patterns[best_j], iupac_patterns[best_i], best_comp,
                                                        bg_model->getV()[0], best_shift);
      }
      else {
        merged_pattern = new IUPACPattern(iupac_patterns[best_i], iupac_patterns[best_j], best_comp,
                                                        bg_model->getV()[0], best_shift);
      }

      std::cout << "merge: " << iupac_patterns[best_j]->get_pattern_string() <<
          " + " << iupac_patterns[best_i]->get_pattern_string() <<
          " -> " << merged_pattern->get_pattern_string() << std::endl;

      delete iupac_patterns[best_j];
      delete iupac_patterns[best_i];

      iupac_patterns.erase(iupac_patterns.begin() + best_j);
      iupac_patterns.erase(iupac_patterns.begin() + best_i);

      iupac_patterns.push_back(merged_pattern);

      //continue searching for merges
      found_better = true;
    }
  }
}

void Peng::process(const int max_pattern_length, const float zscore_threshold, const size_t count_threshold,
                   const int pseudo_counts, const OPTIMIZATION_SCORE opt_score_type, const bool use_em,
                   const float em_saturation_factor, const float min_em_threshold, const int em_max_iterations,
                   const bool use_merging, const float bit_factor_merge_threshold,
                   const bool adv_pwm,
                   std::vector<IUPACPattern*>& best_iupac_patterns) {

  if(max_pattern_length > log(SIZE_MAX) / log(IUPAC_ALPHABET_SIZE) - 1 ||
      max_pattern_length > log(SIZE_MAX) / log(Alphabet::getSize()) - 1) {
    std::cerr << "Warning: pattern length too long!" << std::endl;
    std::cerr << "max pattern length: " << std::max(log(SIZE_MAX) / log(IUPAC_ALPHABET_SIZE) - 1, log(SIZE_MAX) / log(Alphabet::getSize()) - 1) << std::endl;
    exit(1);
  }

  for(size_t pattern_length = std::min(6, max_pattern_length); pattern_length <= max_pattern_length; pattern_length += 2) {
    BasePattern* base_pattern = new BasePattern(pattern_length, strand, k, max_k, sequence_set, bg_model);
    size_t* pattern_counter = base_pattern->getPatternCounter();

    std::vector<size_t> filtered_base_patterns;
    base_pattern->filter_base_patterns(zscore_threshold, count_threshold, filtered_base_patterns);

    std::vector<IUPACPattern*> unoptimized_iupac_patterns;
    optimize_iupac_patterns(opt_score_type, base_pattern, filtered_base_patterns, unoptimized_iupac_patterns);

    filter_iupac_patterns(pattern_length, unoptimized_iupac_patterns);
    for(auto pattern : unoptimized_iupac_patterns) {
      std::cout << "selected iupac pattern: " << IUPACPattern::toString(pattern->get_pattern(), pattern_length) << std::endl;
    }

    #pragma omp parallel for
    for(int i = 0; i < unoptimized_iupac_patterns.size(); i++) {
      IUPACPattern* pattern = unoptimized_iupac_patterns[i];
      pattern->count_sites(pattern_counter);

      if(adv_pwm) {
        std::cout << "adv pwm: ";
        pattern->calculate_adv_pwm(base_pattern, pseudo_counts, pattern_counter, bg_model->getV()[0]);
      }
      else {
        std::cout << "def pwm: ";
        pattern->calculate_pwm(base_pattern, pseudo_counts, pattern_counter, bg_model->getV()[0]);
      }
      std::cout << IUPACPattern::toString(pattern->get_pattern(), pattern_length) << " -> " << pattern->get_pattern_string() << std::endl;
    }

    std::cout << "max_k: " << this->max_k << std::endl;
    for(int background = 0; background <= this->max_k; background++) {
      float* pattern_bg_probs = base_pattern->getBackgroundProb(background);
      std::cout << "k: " << background << std::endl;
      std::vector<IUPACPattern*> optimized_patterns;
      if(use_em) {
        em_optimize_pwms(unoptimized_iupac_patterns, base_pattern,
                         em_saturation_factor, min_em_threshold, em_max_iterations,
                         pattern_bg_probs, optimized_patterns);
      }

      if(use_merging) {
        if(pattern_length >= MIN_MERGE_OVERLAP) {
            merge_iupac_patterns(pattern_length, bit_factor_merge_threshold, bg_model,
                                 optimized_patterns);
        }
        else {
          std::cerr << "Warning: Specified pattern length ("
              << pattern_length << ") is too low for merging!" << std::endl;
        }
      }

      for(int i = 0; i < optimized_patterns.size(); i++) {
        optimized_patterns[i]->set_optimization_bg_model_order(background);
        best_iupac_patterns.push_back(optimized_patterns[i]);
      }
    }

    delete base_pattern;
  }
}

void Peng::optimize_iupac_patterns(OPTIMIZATION_SCORE score_type,
                                   BasePattern* base_patterns,
                                   std::vector<size_t>& selected_base_patterns,
                                   std::vector<IUPACPattern*>& best_iupac_patterns) {
  std::set<size_t> seen;
  std::set<size_t> best;


  const size_t pattern_length = base_patterns->getPatternLength();
  size_t ltot = base_patterns->getLtot();
  size_t* pattern_counter = base_patterns->getPatternCounter();
  float* pattern_bg = base_patterns->getBackgroundProb(this->k);

  float fraction_pseudo_expected_pattern_counts = 0.1;
  size_t pseudo_expected_pattern_counts = sequence_set->getN() * fraction_pseudo_expected_pattern_counts;

  for(auto pattern : selected_base_patterns) {
    size_t iupac_pattern = base_patterns->baseId2IUPACId(pattern);

    bool found_better_mutant = true;
    float best_score = base_patterns->getOptimizationScore(score_type, pattern, pseudo_expected_pattern_counts);

    IUPACPattern* best_mutant = new IUPACPattern(iupac_pattern, pattern_length);
    best_mutant->calculate_log_pvalue(base_patterns,
                                      ltot,
                                      pattern_bg,
                                      pattern_counter);
    best_mutant->count_sites(pattern_counter);

//    bool already_seen = false;

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
          mutated_pattern->calculate_log_pvalue(base_patterns,
                                                ltot,
                                                pattern_bg,
                                                pattern_counter);
          mutated_pattern->count_sites(pattern_counter);

          //add pattern to currently seen
          current_seen.insert(mutated_pattern->get_pattern());

          float curr_score = mutated_pattern->getOptimizationScore(score_type, ltot, pseudo_expected_pattern_counts);
          if(curr_score < best_score) {
            delete best_mutant;
            found_better_mutant = true;
            best_score = curr_score;
            float exp_count = mutated_pattern->get_bg_p() * ltot;
            best_mutant = mutated_pattern;
            std::cout
				<< "\t" << IUPACPattern::toString(best_mutant->get_pattern(), pattern_length)
            	<< "\t" << best_mutant->getLogPval() << "\t" << (int) exp_count
				<< "\t" << best_mutant->get_sites()
            	<< std::endl;
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

      std::cout << "optimization: " << base_patterns->toString(pattern) << " -> " << IUPACPattern::toString(best_mutant->get_pattern(), pattern_length) << std::endl;
    }
    else {
      std::cout << "optimization: " << base_patterns->toString(pattern) << " removed" << std::endl;
      delete best_mutant;
      best_mutant = nullptr;
    }
  }

  for(auto pattern : best_iupac_patterns) {
    std::cout << "iupac pattern: " << IUPACPattern::toString(pattern->get_pattern(), pattern_length) <<
        "\t" << pattern->get_sites() << "\t" << pattern->getLogPval() << std::endl;
  }
}

void Peng::filter_iupac_patterns(const size_t pattern_length, std::vector<IUPACPattern*>& iupac_patterns) {
  std::vector<IUPACPattern*> deselected_patterns;
  std::vector<IUPACPattern*> kept_patterns;

  for(int i = 0; i < iupac_patterns.size(); i++) {
    IUPACPattern* pat = iupac_patterns[i];
    size_t pattern = pat->get_pattern();

    //count non-informative positions ('N')
    int less_informative_positions = 0;
    for(int p = 0; p < pattern_length; p++) {
      int c = IUPACPattern::getNucleotideAtPos(pattern, p);
      if(c == to_underlying(IUPAC_Alphabet::N)) {
        less_informative_positions += 1;
      }
    }

    //limit fraction of non-informative positions ('N')
    if(pattern_length - less_informative_positions <= 3) {
      deselected_patterns.push_back(pat);
    }
    else {
      kept_patterns.push_back(pat);
    }
  }

  iupac_patterns.clear();
  iupac_patterns.insert(iupac_patterns.begin(), kept_patterns.begin(), kept_patterns.end());
  kept_patterns.clear();

  std::sort(iupac_patterns.begin(), iupac_patterns.end(), sort_IUPAC_patterns);
  float min_pvalue = -5;
  if(iupac_patterns.size() > 0) {
    min_pvalue = std::min(-5.0, iupac_patterns[0]->getLogPval() * 0.2);
  }

  for(int i = 0; i < iupac_patterns.size(); i++) {
    IUPACPattern* pat = iupac_patterns[i];
    if(pat->getLogPval() < min_pvalue) {
      kept_patterns.push_back(pat);
    }
    else{
      deselected_patterns.push_back(pat);
    }
  }

  iupac_patterns.clear();
  iupac_patterns.insert(iupac_patterns.begin(), kept_patterns.begin(), kept_patterns.end());
  kept_patterns.clear();

  for(int i = 0; i < deselected_patterns.size(); i++) {
    delete deselected_patterns[i];
  }
  deselected_patterns.clear();
}


void Peng::printShortMeme(std::vector<IUPACPattern*>& best_iupac_patterns,
                                   const std::string output_filename,
                                   BackgroundModel* bg_model) {

  std::sort(best_iupac_patterns.begin(), best_iupac_patterns.end(), sort_IUPAC_patterns);

  std::ofstream myfile (output_filename);
  if (myfile.is_open()) {
    myfile << "MEME version 4" << std::endl;
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

    for(auto pattern : best_iupac_patterns) {
      myfile << "MOTIF " << pattern->get_pattern_string() << std::endl;
      myfile << "letter-probability matrix:" <<
          " alength= " << 4 <<
          " w= " << pattern->get_pattern_length() <<
          " nsites= " << pattern->get_sites() <<
          " bg_prob= " << pattern->get_bg_p() <<
          " opt_bg_order= " << pattern->get_optimization_bg_model_order() <<
          " log(Pval)= " << pattern->getLogPval()<< std::endl;

      float** pwm = pattern->get_pwm();
      for(size_t w = 0; w < pattern->get_pattern_length(); w++) {
        for(size_t a = 0; a < 4; a++) {
          if(a != 0) {
            myfile << " ";
          }
          myfile << std::fixed << std::setprecision(4) << pwm[w][a];
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

  std::sort(best_iupac_patterns.begin(), best_iupac_patterns.end(), sort_IUPAC_patterns);

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
    for(auto pattern : best_iupac_patterns) {
      myfile << "\t\t{" << std::endl;
      myfile << "\t\t\t\"iupac_motif\" : \"" << pattern->get_pattern_string() << "\"," << std::endl;
      myfile << "\t\t\t\"pattern_length\" : " << pattern->get_pattern_length() << "," << std::endl;
      myfile << "\t\t\t\"sites\" : " << pattern->get_sites() << "," << std::endl;
      myfile << "\t\t\t\"log(Pval)\" : " << pattern->getLogPval() << "," << std::endl;
      myfile << "\t\t\t\"bg_prob\" : " << pattern->get_bg_p() << "," << std::endl;
      myfile << "\t\t\t\"opt_bg_order\" : " << pattern->get_optimization_bg_model_order() << "," << std::endl;
      myfile << "\t\t\t\"pwm\" : [" << std::endl;
      float** pwm = pattern->get_pwm();
      for(size_t w = 0; w < pattern->get_pattern_length(); w++) {
        myfile << "\t\t\t\t\t[";
        for(size_t a = 0; a < 4; a++) {
          myfile << std::fixed << std::setprecision(4) << pwm[w][a];
          if(a != 3) {
            myfile << ", ";
          }
          else{
            myfile << "]";
          }
        }
        if(w != pattern->get_pattern_length() - 1) {
          myfile << ", ";
        }
        myfile << std::endl;
      }
      myfile << "\t\t\t\t]" << std::endl;
      myfile << "\t\t}";
      if(pattern != best_iupac_patterns[best_iupac_patterns.size() - 1]) {
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
