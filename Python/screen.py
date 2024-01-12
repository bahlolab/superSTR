#!/usr/bin/env python

import argparse
import functools
import glob
import numpy as np
import pandas as pd
from mpmath import quadgl, power, binomial, mpmathify
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
from multiprocessing.pool import Pool
from os import path
from os.path import basename



def information_score(row, start_threshold=0.75):
    max_rlen = row.rlen
    info_score = 0.0
    start_idx = int(start_threshold * max_rlen)
    for idx in range(start_idx, max_rlen + 1):
        rlab = str(idx)
        if not np.isnan(row[rlab]):
            info_score += row[rlab] * idx
    return info_score


def infill(if_df, manifest_df):
    all_snames = set(manifest_df.index)
    obs_snames = set(if_df["sample"].unique())
    to_add = all_snames - obs_snames
    rlen = if_df.rlen.max()
    add_dict = {}
    for idx in range(1, rlen + 1):
        add_dict[str(idx)] = 0.0
    add_dict["motif"] = if_df.iloc[0].motif
    add_dict["mlen"] = if_df.iloc[0].mlen
    add_dict["rlen"] = rlen
    for add_samp in to_add:
        add_dict["sample"] = add_samp
        add_dict["grp"] = manifest_df.loc[add_samp].grp
        if_df = if_df.append(add_dict, ignore_index=True)
    return if_df


def cdf(p, x, n):
    cumsum = mpmathify(0)
    for i in range(0, x + 1):
        cumsum = cumsum + binomial(n, i) * power(p, i) * power((1 - p),
                                                               (n - i))
    return cumsum


def p_value_test(n_gt, n_trials, n_grp1, n_grp2):
    m_t = binomial(n_grp1 + n_grp2, n_grp1)
    p_val = mpmathify((n_gt + 1) / (n_trials + 1)) - quadgl(
        lambda x: cdf(x, n_gt, n_trials), [0, (0.5 / m_t)])
    return p_val


def score_df(motif, f_path, manifest_df, control_label, n_trials,
             swap_labs=False, min_threshold=3):
    current_data = pd.read_csv(f_path)
    current_data = current_data.fillna(0.0)
    current_data = infill(current_data, manifest_df)
    current_data["info_score"] = current_data.apply(information_score, axis=1)
    if swap_labs:
        current_data.set_index("sample", inplace=True, drop=False)
        current_data.update(manifest_df)
    targ_list = current_data.grp.unique()
    results = []
    for targ in targ_list:
        if targ == control_label:
            continue
        if len(current_data[current_data["grp"] == targ]) <= min_threshold:
            continue
        cases = current_data.loc[current_data["grp"] == targ]["info_score"]
        controls = current_data.loc[current_data["grp"] != targ]["info_score"]
        if len(cases) == 0 or len(controls) == 0 or sum(cases) == 0.0 or sum(
                controls) == 0.0:
            continue
        mann_whitney_stat, mann_whitney_p = mannwhitneyu(cases, controls,
                                                         alternative="greater")
        n_cases = len(cases)
        n_controls = len(current_data) - n_cases
        mwu_perm_list = list()
        for i in range(0, n_trials):
            current_data = current_data.sample(frac=1.0)
            cases = current_data.iloc[0:n_cases]["info_score"]
            controls = current_data.iloc[n_cases:]["info_score"]
            mwu_result = mannwhitneyu(cases, controls, alternative="greater")
            mwu_perm_list.append(mwu_result[0])
        mwu_perm_list = np.asarray(mwu_perm_list)
        n_gt = (mwu_perm_list > mann_whitney_stat).sum()
        # P-value calculated per Smyth and Phipson
        perm_p = p_value_test(n_gt, n_trials, n_cases, n_controls)
        result_tuple = (
            targ, len(cases), motif, mann_whitney_stat, mann_whitney_p, perm_p)
        results.append(result_tuple)
    return results


def process_file(fname, manifest_path, control_label, use_status, n_perm, min_threshold):
    manifest_df = pd.read_csv(manifest_path, sep="\t",
                              names=["sample", "status", "path", "grp"])
    manifest_df.set_index("sample", inplace=True)
    motif = basename(fname).replace(".csv", "")
    motif_results = score_df(motif, fname, manifest_df, control_label, n_perm,
                             swap_labs=use_status, min_threshold=min_threshold)
    return motif_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--manifest", dest="manifest", help="path to manifest file")
    parser.add_argument("-o", "--output", dest="output", help="output prefix")
    parser.add_argument("-i", "--input", dest="input", help="input path")
    parser.add_argument("--label", dest="df_string", default="grp",
                        help="column label to use for group status")
    parser.add_argument("--controllab", dest="ctrl_label", default="Control",
                        help="label used for controls")
    parser.add_argument("--swaplabs", dest="use_status", action="store_true",
                        help="use manifest file labels")
    parser.add_argument("-@", dest="nprocs", action="store", type=int,
                        help="Number of processes (default 6)", default=6)
    parser.add_argument("--min_motif", dest="min_motif", action="store",
                        type=int,
                        help="Minimum motif length (inclusive, default 3)",
                        default=3)
    parser.add_argument("--max_motif", dest="max_motif", action="store",
                        type=int,
                        help="Maximum motif length (inclusive, default 6)",
                        default=6)
    parser.add_argument("-p", dest="n_permutations", type=int,
                        help="Number of permutations (default 10,000)",
                        default=10000)
    parser.add_argument("--min_thresh", dest="min_thresh", type=int,
                        help="Minimum number of samples for a group to be tested",
                        default=3)
    args = parser.parse_args()
    control_label = args.ctrl_label
    replace_dict = None
    res_list = []
    file_glob_list = []
    for idx in range(args.min_motif, args.max_motif + 1):
        for fname in glob.glob(
                path.join(args.input, "motifs", str(idx) + "mers/*.csv")):
            file_glob_list.append(fname)
    if len(file_glob_list) == 0:
        print("No files loaded! Check that there are files in:")
        for idx in range(args.min_motif, args.max_motif + 1):
            print(path.join(args.input, "motifs", str(idx) + "mers/*.csv"))
        import sys
        sys.exit(1)
    with Pool(processes=args.nprocs) as pool:
        for i in tqdm(pool.imap_unordered(
                functools.partial(process_file, manifest_path=args.manifest,
                                  control_label=control_label,
                                  n_perm=args.n_permutations,
                                  use_status=args.use_status, min_threshold=args.min_thresh), file_glob_list),
                total=len(file_glob_list)):
            for res in i:
                targ, n_cases, motif, mann_whitney_stat, mann_whitney_p, permu_p = res
                res_list.append(
                    {"test class": targ, "N": n_cases, "motif": motif,
                     "MW statistic": mann_whitney_stat,
                     "MW p-value": mann_whitney_p,
                     "MW permutation p-value": permu_p})
    pool.close()
    pool.join()
    res_df = pd.DataFrame(res_list)
    res_df["MW FDR p-value"] = \
        multipletests(res_df["MW permutation p-value"], alpha=0.05,
                      method="fdr_bh")[1]
    res_df.to_csv(args.output, index=False, sep="\t")

