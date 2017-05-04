/*
 * gap_mask.h
 *
 *  Created on: May 3, 2017
 *      Author: mmeier
 */


#include <string.h>
#include <vector>
#include <iostream>

#ifndef SRC_GAP_MASK_H_
#define SRC_GAP_MASK_H_

class GapMask {
 public:
  GapMask(GapMask* init);
  GapMask(int mask_length);
  ~GapMask();
  int get_mask_length();
  bool* get_mask();
  void set_mask(const int pos, const bool is_informative);
  bool is_symmetric();
  bool starts_with_informative();
  bool no_sigle_informative();
  std::string toString();
 private:
  int mask_length;
  bool* mask;
};

void get_masks(const int pattern_length, std::vector<GapMask*>& masks);

#endif /* SRC_GAP_MASK_H_ */
