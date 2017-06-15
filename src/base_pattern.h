/*
 * base_pattern.h
 *
 *  Created on: Jan 30, 2017
 *      Author: mmeier
 */

#ifndef SRC_BASE_PATTERN_H_
#define SRC_BASE_PATTERN_H_

#include <cstdlib>
#include <string>


/**
    BasePattern deals with the encoding of Patterns
    The nucleotide encoding is defined in Alphabet

    The patterns are encoded in numerical values:
    (BaMM-)Alphabet: 'other'<->0, A<->1, C<->2, G<->3, T<->4
    (PEnG-)Alphabet: A<->0, C<->1, G<->2, T<->3
    Alphabet Size a = 4
    Pattern: "ATGC" <-> id: 0*a^0 + 3*a^1 + 2*a^2 + 1*a^3
*/
class BasePattern {
 public:
  //inits factor and pattern_length
  static void init(const size_t pattern_length);

  //returns factors
  static size_t* getFactors();

  //returns pattern_length
  static size_t getPatternLength();

  //get string to pattern id
  static std::string toString(const size_t pattern_id);

  //get pattern id of reverse complementary pattern
  static size_t getRevCompId(const size_t pattern_id);

  //get nucleotide id of pattern at pos
  static int getNucleotideAtPos(const size_t pattern, const size_t pos);
 private:
  // position specific factors used for encoding of base patterns
  static size_t* factor;

  // length of pattern
  static size_t pattern_length;
};

#endif /* SRC_BASE_PATTERN_H_ */
