import glob
import os
import sys
from os.path import basename

import numpy as np
import pandas as pd
from arch.bootstrap import IIDBootstrap


def get_quantile(xval, ci_val=0.95):
    return np.quantile(xval, ci_val)

def bootstrap_ci(df, bootstraps=1000, max_len=150, min_len=90,
                 info_score=False, max_length=False, lb=False, control_lab=None, user_ci = 0.95):
    if info_score:
        df["iscore"] = 0.0
        for i, row in df.iterrows():
            total = 0.0
            for idx in range(min_len, max_len+1):
                total += idx*row[str(idx)]
            df.at[i, "iscore"] = total
        if control_lab:
            score_vec = df[df["grp"]==control_lab]["iscore"]
        else:
            score_vec = df["iscore"]
    elif max_length:
        if control_lab:
            score_vec = df[df["grp"] == control_lab][str(max_len)]
        else:
            score_vec = df[str(max_len)]
    else:
        print("Unsupported operation - select information score or max length to test.")
        import sys
        sys.exit(1)
    if np.sum(score_vec) == 0.0:
        ret_str = "No non-zero values in test"
        return ret_str
    bs = IIDBootstrap(score_vec)
    try:
        ci = bs.conf_int(get_quantile, extra_kwargs={"ci_val":user_ci}, reps=bootstraps, method='bca')
        if lb:
            threshold = ci[0][0]
        else:
            threshold = ci[1][0]
        ret_str = str(threshold)+"\t"
        # If a control label is specified, drop this from the dataframe for
        # testing outliers. Otherwise, test outliers across the dataset.
        if control_lab:
            df = df[~df["grp"].isin([control_lab])]
        if info_score:
            df = df[df["iscore"] >= threshold]
        elif max_length:
            df = df[df[str(max_len)] >= threshold]
        else:
            print("Unsupported operation: select a ")
        ret_str += ",".join(df["sample"].to_list())
        ret_str+= "\t"
        cdict = df["grp"].value_counts().to_dict()
        ret_str += ",".join([f'{key}:{value}' for key, value in cdict.items()])
    except ValueError:
        ret_str = "Failure in outlier detection library"
    return ret_str

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    paths_group = parser.add_argument_group(title="Global options.")
    paths_group.add_argument("--input", action="store", dest="input_path",
                             help="Input path; should be a directory containing "
                                  "one or more superSTR runs.")
    paths_group.add_argument("--output", action="store", dest="output_path", default=None,
                             help="Output path for files. If not set, output to stdout.")
    paths_group.add_argument("--manifest", action="store", dest="manifest",
                             default=None, help="Path to manifest file (optional)")
    paths_group.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                             default=False, help="Enable verbose mode.")
    motif_arguments = parser.add_argument_group(title="Motif targeting options:")
    motif_arguments.add_argument("--pathogenic", action="store_true",
                                    dest="pathos",
                                    default=False, help="Perform outlier "
                                                        "detection in only known pathogenic motifs.")
    motif_arguments.add_argument("--max_motif", action="store",
                                    dest="max_motif",
                                    default=None, help="Set the maximum motif "
                                                       "length to run outlier calculations to. ")
    motif_arguments.add_argument("--motifs", action="store", dest="motifs", default=None,
                             help="Path to a file containing motifs of interest; one line per motif.")
    mode_arguments = parser.add_argument_group(title="Outlier detection methods")
    mode_arguments.add_argument("--bootstrapCI", action="store_true", dest="ci_95",
                                    default=False, help="Generate outlier calls against a specified "
                                    "the 95%% CI estimate of a percentile using the ABC bootstrap")
    ci_arguments = parser.add_argument_group(title="Bootstrapped percentile CI options")
    ci_arguments.add_argument("-ci", action="store", dest="CI", default=0.95,
                              help="Percentile for estimation")
    ci_arguments.add_argument("-lb", action="store_true", dest="use_lb",
                                    default=False, help="Use lower-bound of "
                                    "95%% CI estimate (default is to use upper)")
    ci_arguments.add_argument("--rl", action="store",
                                    dest="read_length", default=None,
                                    help="Manually set read length.")
    ci_arguments.add_argument("--ctrllab", action="store", dest="ctrllab", default=None, help="Label for control data")
    metric_arguments = parser.add_argument_group(title="Metrics and metric options")
    metric_arguments.add_argument("-is", "--info_score", action="store_true", dest="iscore",
                                    default=False, help="Test for outliers in size-weighted repeat score")
    metric_arguments.add_argument("-ml", "--max_length", action="store_true", dest="maxlen",
                                    default=False, help="Test for outliers in max-length reads only")
    metric_arguments.add_argument("--min_len", action="store", dest="min_len", default=None, type=int,
                                    help="Minimum length for information score calculation")
    metric_arguments.add_argument("--max_len", action="store", dest="max_len", default=None, type=int,
                                    help="Maximum length for information score calculation")

    args = parser.parse_args()
    # TODO: Argument sanitising.
    directory_root = args.input_path
    output_path = args.output_path
    if args.pathos:
        known_patho_motifs = ["3mers/AGC.csv", "5mers/AAAGT.csv",
                              "12mers/CCCCGCCCCGCG.csv", "6mers/CCCCGG.csv",
                              "3mers/AAG.csv", "5mers/AATAG.csv",
                              "3mers/ACG.csv", "3mers/CCG.csv",
                              "6mers/AGGCCC.csv", "4mers/AGGC.csv",
                              "5mers/AATGG.csv"]
        filelist = [directory_root+"/motifs/"+x for x in known_patho_motifs]
    elif args.max_motif:
        filelist = []
        for i in range(1, int(args.max_motif)+1):
            dirpath = directory_root+"/motifs/"+str(i)+"mers/"
            if os.path.exists(dirpath) and os.path.isdir(dirpath):
                for x in glob.glob(dirpath+"*.csv"):
                    filelist.append(x)
    else:
        filelist = glob.glob(directory_root+"/motifs/*mers/*.csv")
    filelist = sorted(filelist, key=lambda x: len(x))
    if args.output_path is None:
        f = sys.stdout
    else:
        f = open(args.output_path, 'wt')
    f.write("Motif\tThreshold\tOutlier samples\tGroup counts\n")
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
        if args.ci_95:
            write_string = bootstrap_ci(df, info_score=args.iscore,
                                        max_len=args.max_len, min_len=args.min_len,
                                        lb=args.use_lb, control_lab=args.ctrllab,
                                        user_ci=float(args.CI))
            if write_string is None or "Failure in" in write_string:
                fails.append(basename(file).replace(".csv", ""))
            elif "non-zero" in write_string:
                nonzeros.append(basename(file).replace(".csv", ""))
            else:
                f.write(basename(file).replace(".csv", "") + "\t")
                f.write(write_string)
                f.write("\n")
    if args.verbose:
        if len(fails) != 0:
            print("Failures in outlier detection library (usually due to low variance):", file=sys.stderr)
            for fname in fails:
                print(fname, file=sys.stderr, end=",")
            print("", file=sys.stderr)
        if len(nonzeros) != 0:
            print("Motifs with all zero values in the tested statistics:", file=sys.stderr)
            for fname in nonzeros:
                print(fname, file=sys.stderr, end=",")
            print("", file=sys.stderr)
