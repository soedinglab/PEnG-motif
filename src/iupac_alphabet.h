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
enum class IUPAC_Alphabet { A = 0, C = 1, G = 2, T = 3, S = 4, W = 5, R = 6, Y = 7, M = 8, K = 9, N = 10};

class IUPACAlphabet {
 private:
  static char* base_2_char;
  static int* char_2_base;

  static std::map<int, std::vector<int>> similar_iupac_nucleotides;
  static std::map<int, std::vector<int>> representative_iupac_nucleotides;

 public:
  static void init(char* alphabet);
  static std::vector<int> get_similar_iupac_nucleotides(int c);
  static std::vector<int> get_representative_iupac_nucleotides(int c);
  static char getBase(int c);
  static int getCode(char c);
  static size_t getAlphabetSize();
};



#endif
