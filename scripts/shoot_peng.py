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


def check_executable_presence(executable_name):
    if not shutil.which(executable_name):
        print('|ERROR| Cannot find %s. Please install it and check your PATH variable.'
              % executable_name, file=sys.stderr)
        return False
    return True


RSCRIPT = "evaluateBaMM.R"
PENG = "peng_motif"
FDR = "FDR"


def main():
    parser = argparse.ArgumentParser(description='A wrapper for PEnG that reranks the found motifs')
    parser.add_argument(metavar='FASTA_FILE', dest='fasta_file', type=str,
                        help='file with the input fasta sequences')
    parser.add_argument("-o", metavar='FILE', dest='meme_output_file', type=str,
                        help='best IUPAC motives will be written in FILE in minimal MEME format')
    parser.add_argument("-j", metavar='FILE', dest='json_output_file', type=str,
                        help='best IUPAC motives will be written in OUTPUT_FILE in JSON format')
    parser.add_argument("-d", "--output_directory", metavar='DIR', dest='output_directory', type=str,
                        help='directory for the temporary files')
    parser.add_argument('--background-sequences', metavar='FASTA_FILE', dest='background_sequences', type=str,
                        help='file with fasta sequences to be used for the background model calculation')
    parser.add_argument('-w', metavar='INT', dest='pattern_length', type=int, default=10,
                        help='initial/minimal length of pattern to be searched')
    parser.add_argument('-t', metavar='FLOAT', dest='zscore_threshold', type=float, default=10,
                        help='lower zscore threshold for basic patterns')
    parser.add_argument('--count-threshold', metavar='INT', dest='count_threshold', type=int, default=1,
                        help='lower threshold for counts of base patterns')
    parser.add_argument('--bg-model-order', metavar='INT', dest='bg_model_order', type=int, default=2,
                        help='order of the background model')
    parser.add_argument('--strand', metavar='PLUS|BOTH', dest='strand', type=str, default='BOTH', choices=['PLUS', 'BOTH'],
                        help='select the strand to work on')
    parser.add_argument('--iupac_optimization_score', metavar='LOGPVAL|EXPCOUNTS|MUTUAL_INFO',
                        dest='iupac_optimization_score', type=str, default='LOGPVAL',
                        choices=['ENRICHMENT', 'LOGPVAL', 'MUTUAL_INFO'],
                        help='select iupac optimization score')
    parser.add_argument('--enrich_pseudocount_factor', type=float, default=0.005, metavar="FLOAT",
                        help="add (enrich_pseudocount_factor x #seqs) pseudo counts "
                        "in the EXPCOUNTS optimization")
    parser.add_argument('--no-em', dest='use_em', action='store_false', default=True,
                        help='shuts off the em optimization')
    parser.add_argument('-a', metavar='FLOAT', dest='em_saturation_threshold', type=float, default=1E4,
                        help='saturation factor for em optimization')
    parser.add_argument('--em-threshold', metavar='FLOAT', dest='em_threshold', type=float, default=0.08,
                        help='threshold for finishing the em optimization')
    parser.add_argument('--em-max-iterations', metavar='INT', dest='em_max_iterations', type=int, default=100,
                        help='max number of em optimization iterations')
    parser.add_argument('--no-merging', dest='use_merging', action='store_false', default=True,
                        help='shuts off the merging of patterns')
    parser.add_argument('-b', metavar='FLOAT', dest='bit_factor_threshold', type=float, default=0.4,
                        help='bit factor threshold for merging IUPAC patterns')
    parser.add_argument('--use-default-pwm', action='store_true', dest='use_default_pwm', default=False,
                        help='use the default calculation of the PWM')
    parser.add_argument('--pseudo-counts', metavar='INT', dest='pseudo_counts', type=int, default=10,
                        help='number of pseudo counts for the calculation of the PWM')
    parser.add_argument('--threads', metavar='INT', dest='number_threads', type=float, default=1,
                        help='number of threads to be used for parallelization')
    parser.add_argument('--silent', action='store_true',
                        help='capture and suppress output on stdout')
    parser.add_argument('--no-scoring', action='store_true',
                        help='skip the calculation of the pwm performance score')
    parser.add_argument('--no-neighbor-filtering', action='store_true',
                        help='do not filter similar base patterns before running the optimization')
    parser.add_argument('--minimum-processed-patterns', type=int, default=25,
                        help='minimum number of iupac patterns that are selected for em optimization')

    args = parser.parse_args()

    if args.meme_output_file is None and args.json_output_file is None:
        print("Warning: you did not define an output file (options -o or -j). Stopping here.", file=sys.stderr)
        sys.exit(1)

    required_tools = [PENG]
    if not args.no_scoring:
        required_tools += [RSCRIPT, FDR]

    ready = True
    for tool in required_tools:
        if not check_executable_presence(tool):
            ready = False
    if not ready:
        sys.exit(10)

    output_directory = args.output_directory
    if args.output_directory is None:
        # work in tmp directory
        with tempfile.TemporaryDirectory() as output_directory:
            run_peng(args, output_directory, not args.no_scoring)
    else:
        if not os.path.exists(output_directory):
            # make user defined directory if not exist
            os.makedirs(output_directory)
        run_peng(args, output_directory, not args.no_scoring)


def build_peng_command(args, protected_fasta_file, peng_output_file, peng_json_file):
    command = [PENG, os.path.abspath(protected_fasta_file),
               "-j", os.path.abspath(peng_json_file),
               "-o", os.path.abspath(peng_output_file)]
    if args.background_sequences:
        command += ["--background-sequences", os.path.abspath(args.background_sequences)]
    command += ["-w", str(args.pattern_length)]
    command += ["-t", str(args.zscore_threshold)]
    command += ["--count-threshold", str(args.count_threshold)]
    command += ["--bg-model-order", str(args.bg_model_order)]
    command += ["--strand", args.strand]
    command += ["--iupac_optimization_score", str(args.iupac_optimization_score)]
    command += ["--enrich_pseudocount_factor", str(args.enrich_pseudocount_factor)]
    if not args.use_em:
        command += ["--no-em"]
    command += ["-a", str(args.em_saturation_threshold)]
    command += ["--em-threshold", str(args.em_threshold)]
    command += ["--em-max-iterations", str(args.em_max_iterations)]
    if not args.use_merging:
        command += ["--no-merging"]
    if args.use_default_pwm:
        command += ["--use-default-pwm"]
    command += ["-b", str(args.bit_factor_threshold)]
    command += ["--pseudo-counts", str(args.pseudo_counts)]
    command += ["--threads", str(args.number_threads)]
    command += ['--minimum-processed-patterns', args.minimum_processed_patterns]
    if args.no_neighbor_filtering:
        command.append('--no-neighbor-filtering')

    print(" ".join(str(c) for c in command))
    return [str(c) for c in command]


# FDR -m 10 -k 0 --cvFold 1
def build_fdr_command(args, protected_fasta_file, peng_output_file, output_directory):
    command = [FDR, output_directory, os.path.abspath(protected_fasta_file),
               "--PWMFile", os.path.abspath(peng_output_file)]
    if args.strand == 'PLUS':
        command += ["--ss"]
    command += ["-m", str(10)]
    command += ["-k", str(0)]
    command += ["--cvFold", str(1)]

    print(" ".join(command))
    return command


def run_peng(args, output_directory, run_scoring):
    # FDR takes the prefix from the input fasta-file
    filename, extension = os.path.splitext(args.fasta_file)
    prefix = os.path.basename(filename)
    whitespace_matcher = re.compile(r'\s+')
    prefix = re.sub(whitespace_matcher, '_', prefix)

    peng_output_file = os.path.join(output_directory, prefix + ".tmp.out")
    peng_json_file = os.path.join(output_directory, prefix + ".tmp.json")

    if args.silent:
        stdout = subprocess.DEVNULL
    else:
        stdout = None

    # run peng
    peng_command_line = build_peng_command(args, args.fasta_file, peng_output_file, peng_json_file)
    subprocess.run(peng_command_line, check=True, stdout=stdout)

    with open(peng_json_file) as fh:
        peng_data = json.load(fh)

    if run_scoring:
        # run FDR
        fdr_command_line = build_fdr_command(args, args.fasta_file, peng_output_file, output_directory)
        subprocess.run(fdr_command_line, check=True, stdout=stdout)
        r_output_file = os.path.join(output_directory, prefix + ".bmscore")
        subprocess.run([RSCRIPT, os.path.abspath(output_directory), prefix], check=True, stdout=stdout)

        # run R script
        zoops_scores = dict()
        with open(r_output_file) as fh:
            for line in fh:
                if line.startswith("prefix"):
                    continue
                try:
                    prefix, motif_number, zoops_rank_score, auc5_score, auprc_score, *_ = line.split()
                    motif_number = int(motif_number)
                except ValueError:
                    # either header or weird line. skipping.
                    continue

                try:
                    # note here ausfc score is used for reranking, instead of fract_occ
                    zoops_scores[motif_number] = float(zoops_rank_score)
                except ValueError:
                    zoops_scores[motif_number] = np.nan

        # update information
        patterns = peng_data["patterns"]
        for idx, p in enumerate(patterns):
            if idx + 1 in zoops_scores:
                p["zoops_score"] = zoops_scores[idx + 1]
                print("{} {}".format(p["iupac_motif"], p["zoops_score"]))
            else:
                p["zoops_score"] = np.nan

        peng_data["patterns"] = sorted(peng_data["patterns"], key=lambda k: k['zoops_score'], reverse=True)
    else:
        patterns = peng_data["patterns"]
        for p in patterns:
            p["zoops_score"] = float('nan')

    if args.meme_output_file:
        write_meme(peng_data, args.meme_output_file)
    if args.json_output_file:
        write_json(peng_data, args.json_output_file)


def write_meme(peng_data, peng_output_file):
    with open(peng_output_file, "w") as fh:
        print("MEME version 4", file=fh)
        print(file=fh)

        alphabet_length = peng_data["alphabet_length"]
        print("ALPHABET= " + peng_data["alphabet"], file=fh)
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
            print(
                ("letter-probability matrix: alength= {} w= {} "
                 "nsites= {} bg_prob= {} opt_bg_order= {} log(Pval)= {} zoops_score= {}")
                .format(
                    alphabet_length, p["pattern_length"], p["sites"],
                    p["bg_prob"], p["opt_bg_order"], p["log(Pval)"],
                    p["zoops_score"]
                ), file=fh
            )
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
