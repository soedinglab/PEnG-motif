#!/usr/bin/env python

"""
Created in Feb 2017

@author: Markus Meier
"""

from collections import defaultdict
from enum import IntEnum
import numpy as np
import sys
import argparse

IUPAC_ALPHABET_SIZE = 11


#an enum to encode the int representation of the iupac nucleotides
class IUPACNucleotide(IntEnum):
    A = 0
    C = 1
    G = 2
    T = 3
    S = 4
    W = 5
    R = 6
    Y = 7
    M = 8
    K = 9
    N = 10


#generates the map for the amiguous iupac nucleotides e.g.: N -> A, C, G, T
def init_representative_map():
    representative_iupac_nucleotides = defaultdict(list)

    representative_iupac_nucleotides[IUPACNucleotide.A].append(IUPACNucleotide.A)
    representative_iupac_nucleotides[IUPACNucleotide.C].append(IUPACNucleotide.C)
    representative_iupac_nucleotides[IUPACNucleotide.G].append(IUPACNucleotide.G)
    representative_iupac_nucleotides[IUPACNucleotide.T].append(IUPACNucleotide.T)

    representative_iupac_nucleotides[IUPACNucleotide.S].append(IUPACNucleotide.C)
    representative_iupac_nucleotides[IUPACNucleotide.S].append(IUPACNucleotide.G)

    representative_iupac_nucleotides[IUPACNucleotide.W].append(IUPACNucleotide.A)
    representative_iupac_nucleotides[IUPACNucleotide.W].append(IUPACNucleotide.T)

    representative_iupac_nucleotides[IUPACNucleotide.R].append(IUPACNucleotide.A)
    representative_iupac_nucleotides[IUPACNucleotide.R].append(IUPACNucleotide.G)

    representative_iupac_nucleotides[IUPACNucleotide.Y].append(IUPACNucleotide.C)
    representative_iupac_nucleotides[IUPACNucleotide.Y].append(IUPACNucleotide.T)

    representative_iupac_nucleotides[IUPACNucleotide.M].append(IUPACNucleotide.A)
    representative_iupac_nucleotides[IUPACNucleotide.M].append(IUPACNucleotide.C)

    representative_iupac_nucleotides[IUPACNucleotide.K].append(IUPACNucleotide.G)
    representative_iupac_nucleotides[IUPACNucleotide.K].append(IUPACNucleotide.T)

    representative_iupac_nucleotides[IUPACNucleotide.N].append(IUPACNucleotide.N)
    representative_iupac_nucleotides[IUPACNucleotide.N].append(IUPACNucleotide.N)
    representative_iupac_nucleotides[IUPACNucleotide.N].append(IUPACNucleotide.N)
    representative_iupac_nucleotides[IUPACNucleotide.N].append(IUPACNucleotide.N)

    return representative_iupac_nucleotides


#generates a map to translate the int representation of the iupac nucleotides to chars
def get_iupac_int2char():
    int2char = dict()

    int2char[IUPACNucleotide.A] = 'A'
    int2char[IUPACNucleotide.C] = 'C'
    int2char[IUPACNucleotide.G] = 'G'
    int2char[IUPACNucleotide.T] = 'T'
    int2char[IUPACNucleotide.S] = 'S'
    int2char[IUPACNucleotide.W] = 'W'
    int2char[IUPACNucleotide.R] = 'R'
    int2char[IUPACNucleotide.Y] = 'Y'
    int2char[IUPACNucleotide.M] = 'M'
    int2char[IUPACNucleotide.K] = 'K'
    int2char[IUPACNucleotide.N] = 'N'

    return int2char


# returns a sample bg model; mayhaps better to read from an external file?
def get_bg_model():
    bg_model = np.zeros(4)

    bg_model[IUPACNucleotide.A] = 0.2
    bg_model[IUPACNucleotide.C] = 0.3
    bg_model[IUPACNucleotide.G] = 0.3
    bg_model[IUPACNucleotide.T] = 0.2

    return bg_model


# init the profiles for the iupac nucleotides with the given bg_model
def init_iupac_profiles(representative_iupac_nucleotides, bg_model, c=0.2, t=0.7):
    iupac_profiles = np.zeros((IUPAC_ALPHABET_SIZE, 4), np.float)

    for iupac_c in range(IUPAC_ALPHABET_SIZE):
        rep = representative_iupac_nucleotides[iupac_c]
        for a in range(4):
            iupac_profiles[iupac_c][a] += c * bg_model[a]

            for r in rep:
                if a == r:
                    iupac_profiles[iupac_c][a] += t

    return iupac_profiles


#calculates the distance between two profiles; based on the Shannon Entropy?
def calculate_d(profile1, profile2):
    d = 0.0
    for a in range(4):
        d += (profile1[a] - profile2[a]) * (np.log2(profile1[a]) - np.log2(profile2[a]))
    return d


#finds for each profile in the pwm the closest iupac profile
def get_iupac_string(pwm, iupac_profiles, int2char):
    res = []

    pattern_length = len(pwm)
    for i in range(pattern_length):
        min_dist = np.inf
        min_iupac = 0

        for m in range(IUPAC_ALPHABET_SIZE):
            dist = calculate_d(pwm[i], iupac_profiles[m])

            if dist < min_dist:
                min_dist = dist
                min_iupac = m

        res.append(int2char[min_iupac])

    return "".join(res)


#read the pwm from an external file
def read_pwm(filename):
    pwm = []
    with open(filename) as fh:
        for line in fh:
            profile = np.zeros(4)
            tokens = line.split()

            if len(tokens) != 4:
                print("ERROR: line does not seem to be part of a valid pwm!!!", file=sys.stderr)
                print("\t{}".format(line), file=sys.stderr)
                exit(1)

            for i, token in enumerate(tokens):
                profile[i] = float(token)

            EPSILON = 0.1
            if np.sum(profile) >= 1.0 + EPSILON or np.sum(profile) <= 1.0 - EPSILON:
                print("ERROR: line does not seem to be part of a valid pwm!!!", file=sys.stderr)
                print("\t{}".format(line), file=sys.stderr)
                exit(1)

            pwm.append(profile)
    return pwm


#THE main ;)
def main():
    parser = argparse.ArgumentParser(description='Translates a PWM into an IUPAC identifier and prints it')
    parser.add_argument(metavar='PWM_FILE', dest='pwm_file', type=str,
                        help='file with the pwm')

    args = parser.parse_args()

    #preparation
    representative_iupac_nucleotides = init_representative_map()
    int2char = get_iupac_int2char()
    bg_model = get_bg_model()
    iupac_profiles = init_iupac_profiles(representative_iupac_nucleotides, bg_model)

    #read pwm
    pwm = read_pwm(args.pwm_file)

    #actual translation
    result = get_iupac_string(pwm, iupac_profiles, int2char)
    print(result)


#if called as a script; calls the main method
if __name__ == '__main__':
    main()
