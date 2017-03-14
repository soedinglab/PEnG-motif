#!/usr/bin/env python

"""
Created in Feb 2017

@author: Markus Meier
"""

import argparse
import os
import subprocess
import tempfile
import json
import glob
import numpy as np

RSCRIPT="/home/mmeier/git/bamm-private/R/plotAUSFC_rank.R"
PENG="/home/mmeier/opt/PEnG/bin/peng_motif"
BAMM="/home/mmeier/opt/bamm/BaMMmotif"

def main():
    parser = argparse.ArgumentParser(description='Translates a PWM into an IUPAC identifier and prints it')
    parser.add_argument(metavar='FASTA_FILE', dest='fasta_file', type=str,
                        help='file with the input fasta sequences')
    parser.add_argument("-d", "--output_directory", metavar='DIR', dest='output_directory', type=str,
                        help='file with the input fasta sequences')
    parser.add_argument("-b", "--prefix", metavar='PREFIX', dest='prefix', type=str,
                        help='prefix of output files')

    args = parser.parse_args()

    output_directory = args.output_directory
    if output_directory == None:
        #make temp directory
        output_directory = tempfile.mkdtemp()
    elif not os.path.exists(output_directory):
        #make user defined directory if not exist
        os.makedirs(output_directory)

    prefix = args.prefix
    if prefix == None:
        prefix = ""

    #run peng
    peng_output_file = os.path.join(output_directory, prefix + ".tmp.out")
    peng_json_file = os.path.join(output_directory, prefix + ".tmp.json")
    subprocess.check_output([PENG, args.fasta_file, "-j", peng_json_file,
                                "-o", peng_output_file, "-b", str(10), "-w", str(10), "-a", str(1E5)])

    #run bamm
    subprocess.check_output([BAMM, output_directory, args.fasta_file,
                                "--PWMFile", peng_output_file, "--FDR", "--savePvalues"])

    #run R script
    #benchmark_file = os.path.join(output_directory, prefix + ".tmp.score")
    #./FDRaveragedRecall.R <inputDir> <basename of your file>
    zoops_scores = dict()
    mops_scores = dict()
    for bamm_output in glob.glob(os.path.join(output_directory, prefix + "_motif_*.zoops.pvalues")):
        bn = os.path.basename(bamm_output)
        bn = bn.replace(".zoops.pvalues", "")
        index = int(bn.split("_motif_")[1]) - 1
        print("index {}".format(index))

        try:
            subprocess.check_output(["Rscript", "--vanilla", RSCRIPT, output_directory, bn])
        except:
            continue

        with open(os.path.join(output_directory, bn + ".rankscore")) as fh:
            lines = fh.readlines()

        for idx, line in enumerate(lines):
            if line.startswith("Ranking score (ZOOPS):"):
                if lines[idx+1] != "NA\n":
                    zoops_scores[index] = float(lines[idx+1])
                else:
                    zoops_scores[index] = np.nan
                print("\tzoops {}".format(float(lines[idx+1])))
            elif line.startswith("Ranking score (MOPS):"):
                if lines[idx+1] != "NA\n":
                    mops_scores[index] = float(lines[idx+1])
                else:
                    mops_scores[index] = np.nan
                print("\tmops {}".format(mops_scores[index]))

    with open(peng_json_file) as fh:
        peng_data = json.load(fh)

    #update information
    patterns = peng_data["patterns"]
    for idx, p in enumerate(patterns):
        if idx in zoops_scores:
            p["zoops_score"] = zoops_scores[idx]
        else:
            p["zoops_score"] = np.nan
        if idx in mops_scores:
            p["mops_score"] = mops_scores[idx]
        else:
            p["mops_score"] = np.nan

    peng_data["patterns"] = sorted(peng_data["patterns"], key=lambda k: k['zoops_score'], reverse=True)

    peng_final_output_file = os.path.join(output_directory, prefix + ".out")
    write_meme(peng_data, peng_final_output_file)

    peng_final_json_file = os.path.join(output_directory, prefix + ".json")
    write_json(peng_data, peng_final_json_file)

def write_meme(peng_data, peng_output_file):
    with open(peng_output_file, "w") as fh:
        print("MEME version 4", file=fh)
        print(file=fh)

        alphabet_length = peng_data["alphabet_length"]
        print("ALPHABET= "+peng_data["alphabet"], file=fh)
        print(file=fh)

        print("Background letter frequencies", file=fh)

        bg_probs = []
        for idx, nt in enumerate(peng_data["alphabet"]):
            bg_probs.append(nt)
            bg_probs.append(str(peng_data["bg"][idx]))
        print(" ".join(bg_probs), file=fh)
        print(file=fh)

        patterns = peng_data["patterns"]

        for p in patterns:
            print("MOTIF {}".format(p["iupac_motif"]), file=fh)
            print(("letter-probability matrix: alength= {} w= {} "
                    "nsites= {} bg_prob= {} log(Pval)= {} zoops_score= {} mops_score= {}").format(alphabet_length, p["pattern_length"], p["sites"], p["bg_prob"], p["log(Pval)"], p["zoops_score"], p["mops_score"]), file=fh)
            pwm = p["pwm"]

            for line in pwm:
                print(" ".join(['{:.4f}'.format(x) for x in line]), file=fh)
            print(file=fh)

def write_json(peng_data, json_output_file):
    with open(json_output_file, 'w') as fh:
        json.dump(peng_data, fh)


#if called as a script; calls the main method
if __name__ == '__main__':
    main()
