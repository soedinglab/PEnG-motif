/*
 * base_pattern.cpp
 *
 *  Created on: Jan 30, 2017
 *      Author: mmeier
 */

#include "shared/Alphabet.h"
#include <cstdio>
#include <math.h>
#include "base_pattern.h"

//declare static variables
size_t* BasePattern::factor = nullptr;
size_t BasePattern::pattern_length = 0;

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
  size_t tmp_id = pattern_id;
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

