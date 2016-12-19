/*
 * iupac_alphabet.cpp
 *
 *  Created on: Nov 16, 2016
 *      Author: mmeier
 */

#include "iupac_alphabet.h"

std::map<uint8_t, std::vector<uint8_t>> IUPACAlphabet::similar_iupac_nucleotides;
std::map<uint8_t, std::vector<uint8_t>> IUPACAlphabet::representative_iupac_nucleotides;
char* IUPACAlphabet::base_2_char = NULL;

void IUPACAlphabet::init(char* alphabet) {
  similar_iupac_nucleotides[A] = std::vector<uint8_t>();
  similar_iupac_nucleotides[A].push_back(A);
  similar_iupac_nucleotides[A].push_back(W);
  similar_iupac_nucleotides[A].push_back(R);
  similar_iupac_nucleotides[A].push_back(M);
  similar_iupac_nucleotides[A].push_back(N);

  similar_iupac_nucleotides[C] = std::vector<uint8_t>();
  similar_iupac_nucleotides[C].push_back(C);
  similar_iupac_nucleotides[C].push_back(S);
  similar_iupac_nucleotides[C].push_back(Y);
  similar_iupac_nucleotides[C].push_back(M);
  similar_iupac_nucleotides[C].push_back(N);

  similar_iupac_nucleotides[G] = std::vector<uint8_t>();
  similar_iupac_nucleotides[G].push_back(G);
  similar_iupac_nucleotides[G].push_back(S);
  similar_iupac_nucleotides[G].push_back(R);
  similar_iupac_nucleotides[G].push_back(K);
  similar_iupac_nucleotides[G].push_back(N);

  similar_iupac_nucleotides[T] = std::vector<uint8_t>();
  similar_iupac_nucleotides[T].push_back(T);
  similar_iupac_nucleotides[T].push_back(W);
  similar_iupac_nucleotides[T].push_back(Y);
  similar_iupac_nucleotides[T].push_back(K);
  similar_iupac_nucleotides[T].push_back(N);

  similar_iupac_nucleotides[S] = std::vector<uint8_t>();
  similar_iupac_nucleotides[S].push_back(S);
  similar_iupac_nucleotides[S].push_back(C);
  similar_iupac_nucleotides[S].push_back(G);
  similar_iupac_nucleotides[S].push_back(R);
  similar_iupac_nucleotides[S].push_back(Y);
  similar_iupac_nucleotides[S].push_back(M);
  similar_iupac_nucleotides[S].push_back(K);
  similar_iupac_nucleotides[S].push_back(N);

  similar_iupac_nucleotides[W] = std::vector<uint8_t>();
  similar_iupac_nucleotides[W].push_back(W);
  similar_iupac_nucleotides[W].push_back(A);
  similar_iupac_nucleotides[W].push_back(T);
  similar_iupac_nucleotides[W].push_back(R);
  similar_iupac_nucleotides[W].push_back(Y);
  similar_iupac_nucleotides[W].push_back(M);
  similar_iupac_nucleotides[W].push_back(K);
  similar_iupac_nucleotides[W].push_back(N);

  similar_iupac_nucleotides[R] = std::vector<uint8_t>();
  similar_iupac_nucleotides[R].push_back(R);
  similar_iupac_nucleotides[R].push_back(A);
  similar_iupac_nucleotides[R].push_back(G);
  similar_iupac_nucleotides[R].push_back(S);
  similar_iupac_nucleotides[R].push_back(W);
  similar_iupac_nucleotides[R].push_back(M);
  similar_iupac_nucleotides[R].push_back(K);
  similar_iupac_nucleotides[R].push_back(N);

  similar_iupac_nucleotides[Y] = std::vector<uint8_t>();
  similar_iupac_nucleotides[Y].push_back(Y);
  similar_iupac_nucleotides[Y].push_back(C);
  similar_iupac_nucleotides[Y].push_back(T);
  similar_iupac_nucleotides[Y].push_back(S);
  similar_iupac_nucleotides[Y].push_back(W);
  similar_iupac_nucleotides[Y].push_back(M);
  similar_iupac_nucleotides[Y].push_back(K);
  similar_iupac_nucleotides[Y].push_back(N);

  similar_iupac_nucleotides[M] = std::vector<uint8_t>();
  similar_iupac_nucleotides[M].push_back(M);
  similar_iupac_nucleotides[M].push_back(A);
  similar_iupac_nucleotides[M].push_back(C);
  similar_iupac_nucleotides[M].push_back(S);
  similar_iupac_nucleotides[M].push_back(W);
  similar_iupac_nucleotides[M].push_back(R);
  similar_iupac_nucleotides[M].push_back(Y);
  similar_iupac_nucleotides[M].push_back(N);

  similar_iupac_nucleotides[K] = std::vector<uint8_t>();
  similar_iupac_nucleotides[K].push_back(K);
  similar_iupac_nucleotides[K].push_back(G);
  similar_iupac_nucleotides[K].push_back(T);
  similar_iupac_nucleotides[K].push_back(S);
  similar_iupac_nucleotides[K].push_back(W);
  similar_iupac_nucleotides[K].push_back(R);
  similar_iupac_nucleotides[K].push_back(Y);
  similar_iupac_nucleotides[K].push_back(N);

  similar_iupac_nucleotides[N] = std::vector<uint8_t>();
  similar_iupac_nucleotides[N].push_back(N);
  similar_iupac_nucleotides[N].push_back(A);
  similar_iupac_nucleotides[N].push_back(C);
  similar_iupac_nucleotides[N].push_back(G);
  similar_iupac_nucleotides[N].push_back(T);
  similar_iupac_nucleotides[N].push_back(S);
  similar_iupac_nucleotides[N].push_back(W);
  similar_iupac_nucleotides[N].push_back(R);
  similar_iupac_nucleotides[N].push_back(Y);
  similar_iupac_nucleotides[N].push_back(M);
  similar_iupac_nucleotides[N].push_back(K);

  representative_iupac_nucleotides[A] = std::vector<uint8_t>();
  representative_iupac_nucleotides[A].push_back(A);

  representative_iupac_nucleotides[C] = std::vector<uint8_t>();
  representative_iupac_nucleotides[C].push_back(C);

  representative_iupac_nucleotides[G] = std::vector<uint8_t>();
  representative_iupac_nucleotides[G].push_back(G);

  representative_iupac_nucleotides[T] = std::vector<uint8_t>();
  representative_iupac_nucleotides[T].push_back(T);

  representative_iupac_nucleotides[S] = std::vector<uint8_t>();
  representative_iupac_nucleotides[S].push_back(C);
  representative_iupac_nucleotides[S].push_back(G);

  representative_iupac_nucleotides[W] = std::vector<uint8_t>();
  representative_iupac_nucleotides[W].push_back(A);
  representative_iupac_nucleotides[W].push_back(T);

  representative_iupac_nucleotides[R] = std::vector<uint8_t>();
  representative_iupac_nucleotides[R].push_back(A);
  representative_iupac_nucleotides[R].push_back(G);

  representative_iupac_nucleotides[Y] = std::vector<uint8_t>();
  representative_iupac_nucleotides[Y].push_back(C);
  representative_iupac_nucleotides[Y].push_back(T);

  representative_iupac_nucleotides[M] = std::vector<uint8_t>();
  representative_iupac_nucleotides[M].push_back(A);
  representative_iupac_nucleotides[M].push_back(C);

  representative_iupac_nucleotides[K] = std::vector<uint8_t>();
  representative_iupac_nucleotides[K].push_back(G);
  representative_iupac_nucleotides[K].push_back(T);

  representative_iupac_nucleotides[N] = std::vector<uint8_t>();
  representative_iupac_nucleotides[N].push_back(A);
  representative_iupac_nucleotides[N].push_back(C);
  representative_iupac_nucleotides[N].push_back(G);
  representative_iupac_nucleotides[N].push_back(T);

  //TODO non-standard alphabets

  base_2_char = ( char* )calloc( 128, sizeof( char ) );
  base_2_char[A] = 'A';
  base_2_char[C] = 'C';
  base_2_char[G] = 'G';
  base_2_char[T] = 'T';
  base_2_char[S] = 'S';
  base_2_char[W] = 'W';
  base_2_char[R] = 'R';
  base_2_char[Y] = 'Y';
  base_2_char[M] = 'M';
  base_2_char[K] = 'K';
  base_2_char[N] = 'N';
}

std::vector<uint8_t> IUPACAlphabet::get_similar_iupac_nucleotides(uint8_t c) {
  return similar_iupac_nucleotides[c];
}

std::vector<uint8_t> IUPACAlphabet::get_representative_iupac_nucleotides(uint8_t c) {
  return representative_iupac_nucleotides[c];
}

char IUPACAlphabet::getBase(uint8_t c) {
  return base_2_char[c];
}

size_t IUPACAlphabet::getAlphabetSize() {
  return representative_iupac_nucleotides.size();
}

