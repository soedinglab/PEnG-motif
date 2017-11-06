//
// Created by Christian Roth on 16.08.17.
//

#ifndef PENG_MOTIF_UTILS_H
#define PENG_MOTIF_UTILS_H

#include <cmath>

inline float calculate_mutual_information(float pattern_observed, float pattern_expected,
                                   unsigned int n_sequences, float prior) {
  float p = 1 - exp(- pattern_observed / n_sequences);
  float q = 1 - exp(- pattern_expected / n_sequences);
  float p_x_1 = prior * p + (1-prior) * q;
  float p_x_0 = 1 - p_x_1;
  float mut_info = prior * (p * (log(p) - log(p_x_1)) + (1 - p) * (log(1-p) - log(p_x_0)));
  mut_info += (1-prior) * (q * (log(q) - log(p_x_1)) + (1-q) * (log(1-q) - log(p_x_0)));
  return mut_info;
}

inline float calculate_entropy_base2(float p) {
  return - p*log2(p) - (1-p)*log2(1-p);
}

inline float calculate_entropy(float p) {
  return - p*log(p) - (1-p)*log(1-p);
}

inline float calculate_mutual_information_fast(float pattern_observed, float pattern_expected,
                                          unsigned int n_sequences, float prior) {
  float p_obs = 1 - exp(- pattern_observed / n_sequences);
  float p_exp = 1 - exp(- pattern_expected / n_sequences);
  float q = prior;
  float p = p_obs * q + p_exp * (1-q);
  auto H = calculate_entropy;
  return - q * H(p_obs) - (1-q) * H(p_exp) + H(p);
}

namespace Utils {
  inline void no_zero_pwm(float** raw_pwm, unsigned rows, unsigned cols, unsigned precision) {
    float delta = std::pow(10, -static_cast<int>(precision));
    float epsilon = delta / (1 - 4*delta);
    for(int i = 0; i < rows; i++) {
      for(int j = 0; j < cols; j++) {
        raw_pwm[i][j] += epsilon;
      }
    }
    IUPACPattern::normalize_pwm(rows, raw_pwm);
  }
}

#endif //PENG_MOTIF_UTILS_H
