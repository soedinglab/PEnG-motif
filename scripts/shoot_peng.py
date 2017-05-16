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
import sys
import re
import numpy as np
import shutil

RSCRIPT = "plotAUSFC_rank.R"
PENG = "peng_motif"
BAMM = "BaMMmotif"


def main():
    parser = argparse.ArgumentParser(description='A wrapper for PEnG that reranks the found motifs')
    parser.add_argument(metavar='FASTA_FILE', dest='fasta_file', type=str,
                        help='file with the input fasta sequences')
    parser.add_argument("-o", metavar='FILE', dest='meme_output_file', type=str,
                        help='best UIPAC motives will be written in FILE in minimal MEME format')
    parser.add_argument("-j", metavar='FILE', dest='json_output_file', type=str,
                        help='best UIPAC motives will be written in OUTPUT_FILE in JSON format')
    parser.add_argument("-d", "--output_directory", metavar='DIR', dest='output_directory', type=str,
                        help='directory for the temporary files')
    parser.add_argument('--background-sequences', metavar='FASTA_FILE', dest='background_sequences', type=str,
                        help='file with fasta sequences to be used for the background model calculation')
    parser.add_argument('-w', metavar='INT', dest='pattern_length', type=int, default=10,
                        help='initial/minimal length of pattern to be searched')
    parser.add_argument('-t', metavar='FLOAT', dest='zscore_threshold', type=float, default=100,
                        help='lower zscore threshold for basic patterns')
    parser.add_argument('--bg-model-order', metavar='INT', dest='bg_model_order', type=int, default=2,
                        help='order of the background model')
    parser.add_argument('--strand', metavar='PLUS|BOTH', dest='strand', type=str, default='BOTH', choices=['PLUS', 'BOTH'],
                        help='select the strand to work on')
    parser.add_argument('--no-em', dest='use_em', action='store_false', default=True,
                        help='shuts off the em optimization')
    parser.add_argument('-a', metavar='FLOAT', dest='em_saturation_threshold', type=float, default=1E5,
                        help='saturation factor for em optimization')
    parser.add_argument('--em-threshold', metavar='FLOAT', dest='em_threshold', type=float, default=0.08,
                        help='threshold for finishing the em optimization')
    parser.add_argument('--em-max-iterations', metavar='INT', dest='em_max_iterations', type=int, default=100,
                        help='max number of em optimization iterations')
    parser.add_argument('--no-merging', dest='use_merging', action='store_false', default=True,
                        help='shuts off the merging of patterns')
    parser.add_argument('-b', metavar='FLOAT', dest='bit_factor_threshold', type=float, default=0.5,
                        help='bit factor threshold for merging IUPAC patterns')
    parser.add_argument('--threads', metavar='INT', dest='number_threads', type=float, default=1,
                        help='number of threads to be used for parallelization')

    args = parser.parse_args()

    if args.meme_output_file is None and args.json_output_file is None:
        print("Warning: you did not define an output file (options -o or -j). Stopping here.", file=sys.stderr)
        sys.exit(0)

    output_directory = args.output_directory
    if args.output_directory is None:
        # work in tmp directory
        with tempfile.TemporaryDirectory() as output_directory:
            run_peng(args, output_directory)
    else:
        if not os.path.exists(output_directory):
            # make user defined directory if not exist
            os.makedirs(output_directory)
        run_peng(args, output_directory)


def build_peng_command(args, protected_fasta_file, peng_output_file, peng_json_file):
    command = [PENG, os.path.abspath(protected_fasta_file),
               "-j", os.path.abspath(peng_json_file),
               "-o", os.path.abspath(peng_output_file)]
    if args.background_sequences:
        command += ["--background-sequences", os.path.abspath(args.background_sequences)]
    command += ["-w", str(args.pattern_length)]
    command += ["-t", str(args.zscore_threshold)]
    command += ["--bg-model-order", str(args.bg_model_order)]
    command += ["--strand", args.strand]
    if not args.use_em:
        command += ["--no-em"]
    command += ["-a", str(args.em_saturation_threshold)]
    command += ["--em-threshold", str(args.em_threshold)]
    command += ["--em-max-iterations", str(args.em_max_iterations)]
    if not args.use_merging:
        command += ["--no-merging"]
    command += ["-b", str(args.bit_factor_threshold)]
    command += ["--threads", str(args.number_threads)]

    print(" ".join(command))
    return command


def build_bamm_command(args, protected_fasta_file, peng_output_file, output_directory):
    command = [BAMM, output_directory, os.path.abspath(protected_fasta_file),
                "--PWMFile", os.path.abspath(peng_output_file), "--FDR", "--savePvalues"]
    command += ["-K", str(args.bg_model_order)]
    if args.strand == 'PLUS':
        command += ["--ss"]
    command += ["--zoops"]

    print(" ".join(command))
    return command


def run_peng(args, output_directory):
    # bamm takes the prefix from the input fasta-file
    filename, extension = os.path.splitext(args.fasta_file)
    prefix = os.path.basename(filename)
    whitespace_matcher = re.compile(r'\s+')
    prefix = re.sub(whitespace_matcher, '_', prefix)

    protected_fasta_file = os.path.join(output_directory, prefix + ".fasta")
    shutil.copyfile(args.fasta_file, protected_fasta_file)

    peng_output_file = os.path.join(output_directory, prefix + ".tmp.out")
    peng_json_file = os.path.join(output_directory, prefix + ".tmp.json")

    # run peng
    peng_command_line = build_peng_command(args, protected_fasta_file, peng_output_file, peng_json_file)
    peng_ret = subprocess.check_output(peng_command_line)

    # run bamm
    bamm_command_line = build_bamm_command(args, protected_fasta_file, peng_output_file, output_directory)
    subprocess.check_output(bamm_command_line)

    r_output_file = os.path.join(output_directory, prefix + ".rank.out")

    subprocess.check_output(["Rscript", "--vanilla", RSCRIPT, os.path.abspath(output_directory), prefix, os.path.abspath(r_output_file)])

    # run R script
    zoops_scores = dict()
    with open(r_output_file) as fh:
        for line in fh:
            if line.startswith("prefix"):
                continue
            prefix, motif_number, zoops_rank_score, fract_occ = line.split()
            motif_number = int(motif_number)

            try:
                zoops_scores[motif_number] = float(zoops_rank_score)
            except:
                zoops_scores[motif_number] = np.nan

    with open(peng_json_file) as fh:
        peng_data = json.load(fh)

    # update information
    patterns = peng_data["patterns"]
    for idx, p in enumerate(patterns):
        if idx + 1 in zoops_scores:
            p["zoops_score"] = zoops_scores[idx + 1]
        else:
            p["zoops_score"] = np.nan

    peng_data["patterns"] = sorted(peng_data["patterns"], key=lambda k: k['zoops_score'], reverse=True)

    if args.meme_output_file:
        write_meme(peng_data, args.meme_output_file)
    if args.json_output_file:
        write_json(peng_data, args.json_output_file)


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
                   "nsites= {} bg_prob= {} log(Pval)= {} zoops_score= {}").format(
                   alphabet_length, p["pattern_length"], p["sites"], p["bg_prob"], p["log(Pval)"],
                   p["zoops_score"]), file=fh)
            pwm = p["pwm"]

            for line in pwm:
                print(" ".join(['{:.4f}'.format(x) for x in line]), file=fh)
            print(file=fh)


def write_json(peng_data, json_output_file):
    with open(json_output_file, 'w') as fh:
        json.dump(peng_data, fh)


# if called as a script; calls the main method
if __name__ == '__main__':
    main()
