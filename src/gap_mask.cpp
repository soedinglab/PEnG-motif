/*
 * gap_mask.cpp
 *
 *  Created on: May 3, 2017
 *      Author: mmeier
 */

#include "gap_mask.h"

#include <sstream>

GapMask::GapMask(GapMask* init) {
  mask_length = init->mask_length;
  mask = new bool[this->mask_length];

  for(int i = 0; i < mask_length; i++) {
    mask[i] = init->mask[i];
  }
}

GapMask::GapMask(int mask_length) {
  this->mask_length = mask_length;
  this->mask = new bool[this->mask_length];
}

GapMask::~GapMask() {
  delete[] mask;
}

int GapMask::get_mask_length() {
  return mask_length;
}

bool* GapMask::get_mask() {
  return mask;
}

void GapMask::set_mask(const int pos, const bool is_informative){
  mask[pos] = is_informative;
}

bool GapMask::is_symmetric() {
  for(int i = 0; i < mask_length; i++) {
    if(mask[i] != mask[mask_length - 1 - i]) {
      return false;
    }
  }
  return true;
}

bool GapMask::starts_with_informative() {
  return mask[0];
}

bool GapMask::no_sigle_informative() {
  int count_inf = 0;
  for (int i = 0; i < mask_length; i++) {
    if(!mask[i]) {
      if(count_inf == 1) {
        return false;
      }
      count_inf = 0;
    }
    else {
      count_inf++;
    }
  }
  return true;
}

std::string GapMask::toString() {
  std::stringstream out;
  out << mask_length;
  out << "\t";
  for (int i = 0; i < mask_length; i++) {
    out << mask[i];
  }
  return out.str();
}


void permutate_mask(const int left_informative, const int left_gaps, GapMask* curr_mask, std::vector<GapMask*>& masks, std::vector<GapMask*>& rejected_masks) {
  if(left_informative == 0 && left_gaps == 0) {
    if(curr_mask->starts_with_informative() &&
        curr_mask->is_symmetric() &&
        curr_mask->no_sigle_informative()) {
      masks.push_back(curr_mask);
    }
    else {
      rejected_masks.push_back(curr_mask);
    }

    return;
  }

  if(left_informative > 0) {
    GapMask* new_mask = new GapMask(curr_mask);
    new_mask->set_mask(new_mask->get_mask_length() - (left_informative + left_gaps), true);
    permutate_mask(left_informative - 1, left_gaps, new_mask, masks, rejected_masks);
  }
  if(left_gaps > 0) {
    curr_mask->set_mask(curr_mask->get_mask_length() - (left_informative + left_gaps), false);
    permutate_mask(left_informative, left_gaps - 1, curr_mask, masks, rejected_masks);
  }
}


void get_masks(const int pattern_length, std::vector<GapMask*>& masks) {
  std::vector<GapMask*> rejected_masks;
  const int max_gaps = pattern_length / 2 + 1;
  for(int nr_gaps = 0; nr_gaps <= max_gaps; nr_gaps++) {
    GapMask* mask = new GapMask(pattern_length + nr_gaps);
    permutate_mask(pattern_length, nr_gaps, mask, masks, rejected_masks);
  }

  for(int i = 0; i < rejected_masks.size(); i++) {
    delete rejected_masks[i];
  }
  rejected_masks.clear();
}
