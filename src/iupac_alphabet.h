/*
 * iupac_alphabet.h
 *
 *  Created on: Nov 16, 2016
 *      Author: mmeier
 */

#ifndef SRC_IUPAC_ALPHABET_H_
#define SRC_IUPAC_ALPHABET_H_

#include <map>
#include <vector>
#include <cstdint>

const int IUPAC_ALPHABET_SIZE = 11;
enum IUPAC_Alphabet { A = 0, C = 1, G = 2, T = 3, S = 4, W = 5, R = 6, Y = 7, M = 8, K = 9, N = 10};

class IUPACAlphabet {
 private:
  static char* base_2_char;
  static uint8_t* char_2_base;

  static std::map<uint8_t, std::vector<uint8_t>> similar_iupac_nucleotides;
  static std::map<uint8_t, std::vector<uint8_t>> representative_iupac_nucleotides;

 public:
  static void init(char* alphabet);
  static std::vector<uint8_t> get_similar_iupac_nucleotides(uint8_t c);
  static std::vector<uint8_t> get_representative_iupac_nucleotides(uint8_t c);
  static char getBase(uint8_t c);
  static uint8_t getCode(char c);
  static size_t getAlphabetSize();
};



#endif
