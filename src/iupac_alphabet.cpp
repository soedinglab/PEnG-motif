/*
 * iupac_alphabet.cpp
 *
 *  Created on: Nov 16, 2016
 *      Author: mmeier
 */

#include "iupac_alphabet.h"
#include "helper-inl.h"
#include <stdlib.h>

std::map<int, std::vector<int>> IUPACAlphabet::similar_iupac_nucleotides {
  {to_underlying(IUPAC_Alphabet::A),
    {
        to_underlying(IUPAC_Alphabet::W),
        to_underlying(IUPAC_Alphabet::R),
        to_underlying(IUPAC_Alphabet::M),
        to_underlying(IUPAC_Alphabet::N)
    }
  },
  {to_underlying(IUPAC_Alphabet::C),
      {
          to_underlying(IUPAC_Alphabet::S),
          to_underlying(IUPAC_Alphabet::Y),
          to_underlying(IUPAC_Alphabet::M),
          to_underlying(IUPAC_Alphabet::N)
      }
  },
  {to_underlying(IUPAC_Alphabet::G),
      {
          to_underlying(IUPAC_Alphabet::S),
          to_underlying(IUPAC_Alphabet::R),
          to_underlying(IUPAC_Alphabet::K),
          to_underlying(IUPAC_Alphabet::N)
      }
  },
  {to_underlying(IUPAC_Alphabet::T),
      {
          to_underlying(IUPAC_Alphabet::W),
          to_underlying(IUPAC_Alphabet::Y),
          to_underlying(IUPAC_Alphabet::K),
          to_underlying(IUPAC_Alphabet::N)
      }
  },
  {to_underlying(IUPAC_Alphabet::S),
      {
          to_underlying(IUPAC_Alphabet::C),
          to_underlying(IUPAC_Alphabet::G),
          to_underlying(IUPAC_Alphabet::R),
          to_underlying(IUPAC_Alphabet::Y),
          to_underlying(IUPAC_Alphabet::M),
          to_underlying(IUPAC_Alphabet::K),
          to_underlying(IUPAC_Alphabet::N)
      }
  },
  {to_underlying(IUPAC_Alphabet::W),
      {
          to_underlying(IUPAC_Alphabet::A),
          to_underlying(IUPAC_Alphabet::T),
          to_underlying(IUPAC_Alphabet::R),
          to_underlying(IUPAC_Alphabet::Y),
          to_underlying(IUPAC_Alphabet::M),
          to_underlying(IUPAC_Alphabet::K),
          to_underlying(IUPAC_Alphabet::N)
      }
  },
  {to_underlying(IUPAC_Alphabet::R),
      {
          to_underlying(IUPAC_Alphabet::A),
          to_underlying(IUPAC_Alphabet::G),
          to_underlying(IUPAC_Alphabet::S),
          to_underlying(IUPAC_Alphabet::W),
          to_underlying(IUPAC_Alphabet::M),
          to_underlying(IUPAC_Alphabet::K),
          to_underlying(IUPAC_Alphabet::N)
      }
  },
  {to_underlying(IUPAC_Alphabet::Y),
      {
          to_underlying(IUPAC_Alphabet::C),
          to_underlying(IUPAC_Alphabet::T),
          to_underlying(IUPAC_Alphabet::S),
          to_underlying(IUPAC_Alphabet::W),
          to_underlying(IUPAC_Alphabet::M),
          to_underlying(IUPAC_Alphabet::K),
          to_underlying(IUPAC_Alphabet::N)
      }
  },
  {to_underlying(IUPAC_Alphabet::M),
      {
          to_underlying(IUPAC_Alphabet::A),
          to_underlying(IUPAC_Alphabet::C),
          to_underlying(IUPAC_Alphabet::S),
          to_underlying(IUPAC_Alphabet::W),
          to_underlying(IUPAC_Alphabet::R),
          to_underlying(IUPAC_Alphabet::Y),
          to_underlying(IUPAC_Alphabet::N)
      }
  },
  {to_underlying(IUPAC_Alphabet::K),
      {
          to_underlying(IUPAC_Alphabet::G),
          to_underlying(IUPAC_Alphabet::T),
          to_underlying(IUPAC_Alphabet::S),
          to_underlying(IUPAC_Alphabet::W),
          to_underlying(IUPAC_Alphabet::R),
          to_underlying(IUPAC_Alphabet::Y),
          to_underlying(IUPAC_Alphabet::N)
      }
  },
  {to_underlying(IUPAC_Alphabet::N),
      {
          to_underlying(IUPAC_Alphabet::A),
          to_underlying(IUPAC_Alphabet::C),
          to_underlying(IUPAC_Alphabet::G),
          to_underlying(IUPAC_Alphabet::T),
          to_underlying(IUPAC_Alphabet::S),
          to_underlying(IUPAC_Alphabet::W),
          to_underlying(IUPAC_Alphabet::R),
          to_underlying(IUPAC_Alphabet::Y),
          to_underlying(IUPAC_Alphabet::M),
          to_underlying(IUPAC_Alphabet::K),
      }
  }
};

std::map<int, std::vector<int>> IUPACAlphabet::representative_iupac_nucleotides {
  {to_underlying(IUPAC_Alphabet::A),
    {
        to_underlying(IUPAC_Alphabet::A)
    }
  },
  {to_underlying(IUPAC_Alphabet::C),
      {
          to_underlying(IUPAC_Alphabet::C)
      }
  },
  {to_underlying(IUPAC_Alphabet::G),
      {
          to_underlying(IUPAC_Alphabet::G)
      }
  },
  {to_underlying(IUPAC_Alphabet::T),
      {
          to_underlying(IUPAC_Alphabet::T)
      }
  },
  {to_underlying(IUPAC_Alphabet::S),
      {
          to_underlying(IUPAC_Alphabet::C),
          to_underlying(IUPAC_Alphabet::G)
      }
  },
  {to_underlying(IUPAC_Alphabet::W),
      {
          to_underlying(IUPAC_Alphabet::A),
          to_underlying(IUPAC_Alphabet::T)
      }
  },
  {to_underlying(IUPAC_Alphabet::R),
      {
          to_underlying(IUPAC_Alphabet::A),
          to_underlying(IUPAC_Alphabet::G)
      }
  },
  {to_underlying(IUPAC_Alphabet::Y),
      {
          to_underlying(IUPAC_Alphabet::C),
          to_underlying(IUPAC_Alphabet::T)
      }
  },
  {to_underlying(IUPAC_Alphabet::M),
      {
          to_underlying(IUPAC_Alphabet::A),
          to_underlying(IUPAC_Alphabet::C),
      }
  },
  {to_underlying(IUPAC_Alphabet::K),
      {
          to_underlying(IUPAC_Alphabet::G),
          to_underlying(IUPAC_Alphabet::T),
      }
  },
  {to_underlying(IUPAC_Alphabet::N),
      {
          to_underlying(IUPAC_Alphabet::A),
          to_underlying(IUPAC_Alphabet::C),
          to_underlying(IUPAC_Alphabet::G),
          to_underlying(IUPAC_Alphabet::T),
      }
  }
};

char* IUPACAlphabet::base_2_char;
int* IUPACAlphabet::char_2_base;

void IUPACAlphabet::init(char* alphabet) {


  base_2_char = ( char* )calloc( 128, sizeof( char ) );
  base_2_char[to_underlying(IUPAC_Alphabet::A)] = 'A';
  base_2_char[to_underlying(IUPAC_Alphabet::C)] = 'C';
  base_2_char[to_underlying(IUPAC_Alphabet::G)] = 'G';
  base_2_char[to_underlying(IUPAC_Alphabet::T)] = 'T';
  base_2_char[to_underlying(IUPAC_Alphabet::S)] = 'S';
  base_2_char[to_underlying(IUPAC_Alphabet::W)] = 'W';
  base_2_char[to_underlying(IUPAC_Alphabet::R)] = 'R';
  base_2_char[to_underlying(IUPAC_Alphabet::Y)] = 'Y';
  base_2_char[to_underlying(IUPAC_Alphabet::M)] = 'M';
  base_2_char[to_underlying(IUPAC_Alphabet::K)] = 'K';
  base_2_char[to_underlying(IUPAC_Alphabet::N)] = 'N';

  char_2_base = ( int* )calloc( 128, sizeof( int ) );
  char_2_base['A'] = to_underlying(IUPAC_Alphabet::A);
  char_2_base['C'] = to_underlying(IUPAC_Alphabet::C);
  char_2_base['G'] = to_underlying(IUPAC_Alphabet::G);
  char_2_base['T'] = to_underlying(IUPAC_Alphabet::T);
  char_2_base['S'] = to_underlying(IUPAC_Alphabet::S);
  char_2_base['W'] = to_underlying(IUPAC_Alphabet::W);
  char_2_base['R'] = to_underlying(IUPAC_Alphabet::R);
  char_2_base['Y'] = to_underlying(IUPAC_Alphabet::Y);
  char_2_base['M'] = to_underlying(IUPAC_Alphabet::M);
  char_2_base['K'] = to_underlying(IUPAC_Alphabet::K);
  char_2_base['N'] = to_underlying(IUPAC_Alphabet::N);
}

std::vector<int> IUPACAlphabet::get_similar_iupac_nucleotides(int c) {
  return similar_iupac_nucleotides[c];
}

std::vector<int> IUPACAlphabet::get_representative_iupac_nucleotides(int c) {
  return representative_iupac_nucleotides[c];
}

char IUPACAlphabet::getBase(int c) {
  return base_2_char[c];
}

int IUPACAlphabet::getCode(char c) {
  return char_2_base[c];
}

size_t IUPACAlphabet::getAlphabetSize() {
  return representative_iupac_nucleotides.size();
}

