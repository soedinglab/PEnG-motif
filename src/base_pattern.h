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

    The pattern encoding is a "binary" representation:
    Alphabet: 'other'<->0, A<->1, C<->2, G<->3, T<->4
    We do not need 'other': A<->0, C<->1, G<->2, T<->3
    Alphabet Size a = 4
    Pattern: ATGC <-> 0*a^0 + 3*a^1 + 2*a^2 + 1*a^3
*/
class BasePattern {
 public:
  /**
      Inits static variables of BasePattern
      Needs to be called before BasePattern is used

      @param pattern_length length of base patterns
  */
  static void init(const size_t pattern_length);

  /**
      Getter for factors used to encode base pattern ids
  */
  static size_t* getFactors();

  /**
      Getter for pattern length
  */
  static size_t getPatternLength();

  /**
      Get string of encoded base pattern id
      @param pattern_id encoded base pattern id
      @return string of base pattern
  */
  static std::string toString(const size_t pattern_id);
  static std::string toString(size_t pattern_id, const int pattern_length, size_t* factors);

  /**
      Get id of corresponding pattern on the minus strand
      @param pattern_id encoded base pattern id (on the plus strand)
      @return encoded pattern id on the minus strand
  */
  static size_t getMinusId(const size_t pattern_id);

  /**
      get nucleotide id of pattern id at position pos
      @param pattern encoded base pattern id
      @param pos position in the pattern
      @return encoded nucleotide
  */
  static int getNucleotideAtPos(const size_t pattern, const size_t pos);
  static int getNucleotideAtPos(const size_t pattern, const size_t pos, size_t* factor);
 private:
  // position specific factors used for encoding base patterns
  static size_t* factor;

  // length of pattern
  static size_t pattern_length;
};

#endif /* SRC_BASE_PATTERN_H_ */
