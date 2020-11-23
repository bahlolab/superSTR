import itertools
from math import floor, ceil
from os.path import basename

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.ticker as tkr
import numpy as np
from matplotlib.colors import ListedColormap

def plot_single_info_plot(g, num_classes, write_path, title, hist=False, kde=True,
                          figsz=(6,10), clab="Control", colmap=None, ncol=1, sharex=False):

    formatter = tkr.ScalarFormatter(useMathText=True, useOffset=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    all_totals = []
    share_ax = None
    totals_dict = {}
    minval = 100000000000
    maxval = -10000000000000
    for i, (ind, grp) in enumerate(g):
        if len(grp) == 1:
            continue
        grp_lab = next(iter(grp["grp"]))
        pltgrp = grp.drop("grp", axis=1)
        print(i)
        print(ind)
        print(pltgrp)
        pltgrp.columns = pltgrp.columns.astype(int)
        passtot=0
        count = 0
        totals = []
        zero_vals = 0
        for line in pltgrp.iterrows():
            count += 1
            total = 0.0
            for i in pltgrp.columns:
                if not np.isnan(line[1][i]):
                    total += i*line[1][i]
            if total == 0.0:
                zero_vals += 1
                totals.append(0)
                all_totals.append(0)
            else:
                totals.append(total)
                all_totals.append(total)
        totals_dict[grp_lab] = totals
        minval = min(minval, min(totals))
        maxval = max(maxval, max(totals))
    colnames = ["value", "category", "disease"]
    df = pd.DataFrame(columns=colnames)
    for disease in totals_dict:
        print(disease, len(totals_dict[disease]))
        if disease == CTRL_STR:
            continue
        else:
            for i in totals_dict[CTRL_STR]:
                df = df.append({"value":float(i), "category":"Diversity", "disease":disease}, ignore_index=True)
            for i in totals_dict[disease]:
                df = df.append({"value":float(i), "category":"Affected", "disease":disease}, ignore_index=True)
    i=0
    i += 1
    plot_df = df
    f = plt.figure(constrained_layout=True, figsize=figsz)
    nrows = ceil((1.0*len(colmap)-1)/ncol)
    spec = gridspec.GridSpec(ncols=ncol, nrows=nrows, figure=f)
    idx=0
    shared_ax = None
    for disease in totals_dict:
        if disease == "Control" or disease == CTRL_STR:
            continue
        if idx == 0:
            xidx = 0
            yidx = 0
            ax = f.add_subplot(spec[0, 0])
            if sharex:
                shared_ax = ax
        else:
            xidx = idx%nrows
            yidx = idx//nrows
            ax = f.add_subplot(spec[xidx, yidx], sharex=shared_ax)
        idx+=1
        aff_col = colmap[disease]
        control_col = colmap[CTRL_STR]
        curr_plot = plot_df[plot_df["disease"].isin([disease])]
        vp = sns.violinplot(y="disease", x="value", hue="category", linewidth=0.0,
                        hue_order=["Affected","Diversity"], scale_hue=False,
                        data=curr_plot, ax=ax,
                        palette={"Affected":aff_col, "Diversity":control_col},
                        split=True, orient='h', inner=None, cut=0.1)
    #     vp.set_ylabel(None)
    #     vp.get_legend().remove()
    #     if xidx != nrows and sharex:
    #         plt.setp(vp.get_xticklabels(), visible=False)
    #         plt.setp(vp.get_xaxis(), visible=False)
    #         vp.set_frame_on(False)
    #     else:
    #         vp.set_xlabel("Summary score")
    #     vp.set_yticklabels([disease + "\naffected n = " + str(curr_plot[curr_plot["category"]=="Affected"].shape[0])                            ])
    # sns.despine(top=True, right=True, left=True)
    # if sharex:
    #     shared_ax.xaxis.set_major_formatter(formatter)
    #     left, right = shared_ax.get_xlim()
    #     shared_ax.set_xlim([0.0, right])
    plt.savefig(write_path, bbox_inches="tight")
    # plt.show()
    plt.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    paths_group = parser.add_argument_group(title="Global options.")
    paths_group.add_argument("--input", action="store", dest="input_path",
                             help="Input path; should be a directory containing "
                                  "one or more superSTR runs.")
    paths_group.add_argument("--output", action="store", dest="out_path", help="Output path; should be an existing directory")
    paths_group.add_argument("-rl", action="store", dest="read_length", help="Read length in nt.")
    paths_group.add_argument("--ctrl", action="store", dest="ctrl_str", help="Control string.")
    args = parser.parse_args()
    CTRL_STR = args.ctrl_str
    output_path = args.out_path
    if output_path[len(output_path)-1] != "/":
        output_path = output_path + "/"
    known_patho_motifs = ["3mers/AGC.csv", "5mers/AAAGT.csv",
                          "12mers/CCCCGCCCCGCG.csv", "6mers/CCCCGG.csv",
                          "3mers/AAG.csv", "5mers/AATAG.csv",
                          "3mers/ACG.csv", "3mers/CCG.csv",
                          "6mers/AGGCCC.csv", "4mers/AGGC.csv",
                          "5mers/AATGG.csv"]
    filelist = [args.input_path + "/motifs/" + x for x in known_patho_motifs]
    read_length = int(args.read_length)
    plot_start_len = floor(0.75*read_length)
    plot_end_len = read_length
    for file in filelist:
        motif = basename(file).replace(".csv", "")
        df = pd.read_csv(file, sep=",")
        for i in range(1, read_length + 1):
            if i < plot_start_len or i > plot_end_len:
                df.drop(str(i), axis=1, inplace=True)
        df.drop("mlen", inplace=True, axis=1)
        df.drop("motif", inplace=True, axis=1)
        df.drop("rlen", inplace=True, axis=1)
        df.set_index("sample", inplace=True, drop=True)
        df.fillna(value=0, inplace=True)
        # Manual specification of names vector:
        # names = ['Group 1', 'Group 2', 'Group 3', 'Control']
        names = df["grp"].unique()
        names = sorted(names)
        if "Control" in names:
            names.remove("Control")
            names.append("Control")
        height_dict = {}
        for i in names:
            height_dict[i] = len(df[(df['grp'] == i)])
        for i in height_dict:
            print(i, height_dict[i])
        df['grp'] = pd.Categorical(df['grp'], names, ordered=True)
        for col in df.columns:
            if col != "grp":
                df[col].fillna(value=0, inplace=True)
        g = df.groupby(['grp'])
        pal = sns.diverging_palette(240, 10, as_cmap=True)
        # Manually specify color map:
        # collist = {"Group 1":"#355d4e", "Group 2":"#8ddacc", "Group 3":"#4f9880", 'Control':"lightgrey"}
        # Auto-guess color map:
        pal = itertools.cycle(sns.color_palette("husl", len(names)).as_hex())
        collist = {}
        for name in names:
            collist[name] = next(pal)
        collist[CTRL_STR] = "lightgrey"
        plot_single_info_plot(g, len(names),
                              write_path=output_path+motif+".pdf",
                              title=motif, hist=False, kde=True,
                              figsz=(10, 6), colmap=collist, ncol=1)
