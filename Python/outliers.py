#!/usr/bin/env python

import glob
import os
import sys
from os.path import basename

import numpy as np
import pandas as pd
from arch.bootstrap import IIDBootstrap

def get_quantile(xval, pc_val=0.95):
    return np.quantile(xval, pc_val)

def bootstrap_ci(df, motif, bootstraps=1000, max_len=150, min_len=90,
                 info_score=False, max_length=False, lb=False,
                 control_lab=None, user_pc=0.95, use_failover=False):
    failed_over = False
    if info_score:
        df["iscore"] = 0.0
        for i, row in df.iterrows():
            total = 0.0
            for idx in range(min_len, max_len + 1):
                total += idx * row[str(idx)]
            df.at[i, "iscore"] = total
        if control_lab:
            score_vec = df[df["grp"] == control_lab]["iscore"]
        else:
            score_vec = df["iscore"]
    elif max_length:
        if control_lab:
            score_vec = df[df["grp"] == control_lab][str(max_len)]
        else:
            score_vec = df[str(max_len)]
    else:
        print(
            "Unsupported operation - select information score or max length to test.")
        import sys
        sys.exit(1)
    all_zero = False
    if np.sum(score_vec) == 0.0:
        all_zero = True
        threshold = 0.0
    else:
        bs = IIDBootstrap(score_vec)
        try:
            ci = bs.conf_int(get_quantile, extra_kwargs={"pc_val": user_pc},
                             reps=bootstraps, method='bca')
        except:
            if use_failover:
                print("WARN: error in BCa bootstrap. Failing over to percentile bootstrap on motif " + df["motif"][0])
                ci = bs.conf_int(get_quantile, extra_kwargs={"pc_val": user_pc},
                                 reps=bootstraps, method='percentile')
                failed_over = True
            else:
                ret_str = "Failure"
                return ret_str
        if lb:
            threshold = ci[0][0]
        else:
            threshold = ci[1][0]
    ret_str = str(threshold) + "\t"
    # If a control label is specified, drop this from the dataframe for
    # testing outliers. Otherwise, test outliers across the dataset.
    if control_lab:
        df = df[~df["grp"].isin([control_lab])]
    if info_score:
        df = df[df["iscore"] >= threshold]
    elif max_length:
        df = df[df[str(max_len)] >= threshold]
    else:
        print("Unsupported operation: select a metric for calculation.")
    sample_ids = df["sample"].to_list()
    ret_str += ",".join(sample_ids)
    ret_str += "\t"
    cdict = df["grp"].value_counts().to_dict()
    ret_str += ",".join([f'{key}:{value}' for key, value in cdict.items()])
    if failed_over:
        if use_failover:
            ret_str+="\t%"
        else:
            ret_str!="\tFAILED"
    else:
        if all_zero:
            ret_str+="\tZEROES"
        elif len(sample_ids) == 0:
            ret_str+="\tNONE"
        else:
            ret_str+="\tBCA"
    return ret_str


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    paths_group = parser.add_argument_group(title="Global options.")
    paths_group.add_argument("-i", "--input", action="store", dest="input_path",
                             help="Input path; should be a directory containing "
                                  "one or more superSTR runs.")
    paths_group.add_argument("-o", action="store", dest="output_path",
                             default=None,
                             help="Output path for files. If not set, output to stdout.")
    paths_group.add_argument("-m", action="store", dest="manifest",
                             default=None,
                             help="Path to manifest file (optional)")
    paths_group.add_argument("-v", "--verbose", action="store_true",
                             dest="verbose",
                             default=False, help="Enable verbose mode.")
    motif_arguments = parser.add_argument_group(
        title="Motif targeting options:")
    motif_arguments.add_argument("--pathogenic", action="store_true",
                                 dest="pathos",
                                 default=False, help="Perform outlier "
                                                     "detection in only known pathogenic motifs.")
    motif_arguments.add_argument("--max_motif", action="store",
                                 dest="max_motif", type=int,
                                 default=None, help="Set the maximum motif "
                                                    "length to run outlier calculations to. ")
    motif_arguments.add_argument("--motifs", action="store", dest="motifs",
                                 default=None,
                                 help="Path to a file containing motifs of interest; one line per motif.")
    mode_arguments = parser.add_argument_group(
        title="Outlier detection methods")
    mode_arguments.add_argument("--bootstrapCI", action="store_true",
                                dest="ci_95",
                                default=False,
                                help="Generate outlier calls against a specified "
                                     "the 95%% CI estimate of a percentile using the ABC bootstrap")
    ci_arguments = parser.add_argument_group(
        title="Bootstrapped percentile CI options")
    ci_arguments.add_argument("-pc", action="store", dest="pc", default=0.95,
                              type=float,
                              help="Percentile for estimation; must be between 0.0 and 1.0")
    ci_arguments.add_argument("-lb", action="store_true", dest="use_lb",
                              default=False, help="Use lower-bound of "
                                                  "95%% CI estimate (default is to use upper)")
    ci_arguments.add_argument("--rl", action="store",
                              dest="read_length", default=None, type=int,
                              help="Manually set read length.")
    ci_arguments.add_argument("--controllab", action="store", dest="ctrllab",
                              default=None, help="Label for control data")
    ci_arguments.add_argument("--failover", action="store_true", dest="failover", default=False,
                              help="Automatically fail over to the percentile bootstrap.")
    metric_arguments = parser.add_argument_group(
        title="Metrics and metric options")
    metric_arguments.add_argument("-is", "--info_score", action="store_true",
                                  dest="iscore",
                                  default=False,
                                  help="Test for outliers in size-weighted repeat score")
    metric_arguments.add_argument("-ml", "--max_length", action="store_true",
                                  dest="maxlen",
                                  default=False,
                                  help="Test for outliers in max-length reads only")
    metric_arguments.add_argument("--min_len", action="store", dest="min_len",
                                  default=None, type=int,
                                  help="Minimum length for information score calculation")
    metric_arguments.add_argument("--max_len", action="store", dest="max_len",
                                  default=None, type=int,
                                  help="Maximum length for information score calculation")

    args = parser.parse_args()
    # TODO: Argument sanitising.
    directory_root = args.input_path
    output_path = args.output_path
    if args.pathos:
        known_patho_motifs = ["3mers/AGC.csv", "5mers/AAAGT.csv", "5mers/AAATG.csv",
                              "12mers/CCCCGCCCCGCG.csv", "6mers/CCCCGG.csv",
                              "3mers/AAG.csv", "5mers/AATAG.csv",
                              "3mers/ACG.csv", "3mers/CCG.csv",
                              "6mers/AGGCCC.csv", "4mers/AGGC.csv",
                              "5mers/AATGG.csv"]
        filelist = [directory_root + "/motifs/" + x for x in
                    known_patho_motifs]
    elif args.max_motif:
        filelist = []
        for i in range(1, int(args.max_motif) + 1):
            dirpath = directory_root + "/motifs/" + str(i) + "mers/"
            if os.path.exists(dirpath) and os.path.isdir(dirpath):
                for x in glob.glob(dirpath + "*.csv"):
                    filelist.append(x)
    elif args.motifs:
        filelist = []
        with open(args.motifs, 'rt') as input_motif_file:
            for line in input_motif_file:
                line = line.strip()
                glob_string = directory_root + "/motifs/" + str(
                    len(line)) + "mers/" + line + ".csv"
                if os.path.exists(glob_string):
                    filelist.append(glob_string)
                else:
                    print("WARNING: non-existent path at ", glob_string)
    else:
        filelist = glob.glob(directory_root + "/motifs/*mers/*.csv")
    filelist = sorted(filelist, key=lambda x: len(x))
    if args.output_path is None:
        f = sys.stdout
    else:
        f = open(args.output_path, 'wt')
    f.write("Motif\tThreshold\tOutlier samples\tGroup counts\tStatus\n")
    fails = []
    nonzeros = []
    manifest_dict = None
    if args.manifest:
        manifest_dict = {}
        for line in open(args.manifest, 'rt'):
            spline = line.strip().split("\t")
            manifest_dict[spline[0]] = spline[3]
    for file in filelist:
        df = pd.read_csv(file)
        if args.manifest:
            df.set_index("sample", drop=False, inplace=True)
            df["grp"].update(pd.Series(manifest_dict))
        df.fillna(value=0.0, inplace=True)
        motif_str = basename(file).replace(".csv", "")
        if args.ci_95:
            write_string = bootstrap_ci(df, motif_str, info_score=args.iscore,
                                        max_len=args.max_len, min_len=args.min_len,
                                        lb=args.use_lb, control_lab=args.ctrllab,
                                        user_pc=float(args.pc), use_failover=args.failover)
            if write_string == "Failure":
                fails.append(motif_str)
            else:
                if not args.verbose and not any(v in write_string for v in ("BCA", "%")):
                    continue
                f.write(motif_str + "\t")
                f.write(write_string)
                f.write("\n")
    if len(fails) > 0:
        print("The following motifs experienced failures in outlier detection without failover to percentile bootstraps:")
        for string in fails:
            print(string)

