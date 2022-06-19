#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def information_score(row, start_threshold=0.75, start_idx=None):
    max_rlen = row.rlen
    info_score = 0.0
    start_idx = int(start_threshold*max_rlen)
    for idx in range(start_idx, max_rlen+1):
        rlab = str(idx)
        if not np.isnan(row[rlab]):
            info_score += row[rlab]*idx
    return info_score

def infill(if_df, manifest_df):
    all_snames = set(manifest_df.index)
    obs_snames = set(if_df["sample"].unique())
    to_add = all_snames-obs_snames
    rlen = if_df.rlen.max()
    add_dict = {}
    for idx in range(1, rlen+1):
        add_dict[str(idx)] = 0.0
    add_dict["motif"] = if_df.iloc[0].motif
    add_dict["mlen"] = if_df.iloc[0].mlen
    add_dict["rlen"] = rlen
    for add_samp in to_add:
        add_dict["sample"] = add_samp
        add_dict["grp"] = manifest_df.loc[add_samp].grp
        if_df = if_df.append(add_dict, ignore_index=True)
    return if_df

def import_data(fpath, manifest_path):
    manifest_df = pd.read_csv(manifest_path, sep="\t",
                              names=["sample", "status", "path", "grp"])
    manifest_df.set_index("sample", inplace=True)
    current_data = pd.read_csv(fpath)
    current_data = current_data.fillna(0.0)
    current_data = infill(current_data, manifest_df)
    current_data["info_score"] = current_data.apply(information_score, axis=1)
    current_data.set_index("sample", inplace=True, drop=False)
    current_data.update(manifest_df)
    return current_data

def plot_violin(bg_df, case_df, ax, ctrl_label, case_label, disease_label, aff_col, control_col):
    colnames = ["value", "category", "disease"]
    df = pd.DataFrame(columns=colnames)
    for i, row in bg_df.iterrows():
        df = df.append({"value": row["info_score"], "category": ctrl_label, "disease": disease_label+ "\n("+str(len(case_df))+")"}, ignore_index=True)
    for i, row in case_df.iterrows():
        df = df.append({"value": row["info_score"], "category": case_label ,"disease": disease_label+ "\n("+str(len(case_df))+")"}, ignore_index=True)
    hues = [case_label, ctrl_label]
    pal = {case_label: aff_col, ctrl_label: control_col}
    vp = sns.violinplot(y="disease", x="value", hue="category", linewidth=0.0,
                        hue_order=hues, scale_hue=False,
                        data=df, ax=ax,
                        palette=pal,
                        split=True, orient='h', inner=None, cut=0)
    vp.set_ylabel(None)
    vp.set_xlabel("Information Score")
    vp.get_legend().remove()
    sns.despine(top=True, right=True, left=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", dest="manifest", help="path to manifest file", action="store")
    parser.add_argument("-o", dest="output", help="output prefix; matplotlib will try and save the file in the format specified by the extension (PDF, PNG etc)", action="store")
    parser.add_argument("-i", dest="input", help="input path", action="store")
    parser.add_argument("-fh", dest="fig_height", type=float, action="store",
                        help="figure height (inches) (default: 8.97638)", default=6.49606)
    parser.add_argument("-fw", dest="fig_width",type=float, action="store",
                        help="figure width (inches) (default: 6.49606)", default=8.97638)
    parser.add_argument("--aff_col", dest="aff_col", default="#3B6978", action="store",
                        help="color for affected group; matplotlib color or hex (e.g. 'lightgrey' or '#3B6978'. Defaults to '#3B6978'")
    parser.add_argument("--ctrl_col", dest="ctrl_col", default="lightgrey", action="store",
                        help="color for control/background group; matplotlib color or hex (e.g. 'lightgrey' or '#3B6978'. Defaults to 'lightgrey'")
    parser.add_argument("--ctrl_lab", dest="ctrl_lab", action="store",
                        help="label for control/background group; case sensitive")
    parser.add_argument("--aff_lab", dest="aff_lab", action="store",
                        help="label for affected group; case sensitive")
    parser.add_argument("--downsample_ctrl", dest="downsample_ctrl", default=None, type=float, action="store")
    parser.add_argument("--downsample_case", dest="downsample_case", default=None, type=float, action="store")
    args = parser.parse_args()
    f = plt.figure(figsize=(args.fig_width,args.fig_height))
    ax = plt.gca()
    df = import_data(args.input, args.manifest)
    labs = set(df["grp"].unique())
    if args.ctrl_lab not in labs or args.aff_lab not in labs:
        if args.ctrl_lab not in labs:
            print(args.ctrl_lab, " is not a valid label.")
        if args.aff_lab not in labs:
            print(args.aff_lab, " is not a valid label.")
    bg_df = df[df["grp"] == args.ctrl_lab]
    if args.downsample_ctrl:
        bg_df = bg_df.sample(frag=args.downsample_ctrl)
    aff_df = df[df["grp"] == args.aff_lab]
    if args.downsample_case:
        aff_df = aff_df.sample(frac=args.downsample_case)
    plot_violin(bg_df, aff_df, ax, args.ctrl_lab, args.aff_lab, args.aff_lab, args.aff_col, args.ctrl_col)
    plt.tight_layout()
    plt.savefig(args.output)
