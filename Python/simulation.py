#!/usr/bin/env python

"""Docstring"""
import csv
import gzip
import math
import random
import sys
import zlib
from collections import defaultdict
from pathlib import Path

import numpy
import pandas
import requests
from pyfaidx import Fasta
from requests.adapters import HTTPAdapter
from sklearn import metrics
from sklearn.metrics import auc, precision_recall_curve
from urllib3 import Retry

# DNA alphabet is A/C/G/T
DNA_ALPHABET = ["C", "G", "T", "A"]
# CG proportion of 40.87% per Piovesan et al (2019), BMC Res. Notes)
HSA_DNA_PROBABILITIES = [0.20435, 0.20435, 0.29565, 0.29565]
# AT content of intergenic regions of P falciparum per doi:10.1093/nar/gkw1259
PFA_DNA_PROBABILITIES = [0.068, 0.068, 0.432, 0.432]
# GC content of S. coelicolor per https://doi.org/10.1038/417141a
SCO_DNA_PROBABILITIES = [0.3606, 0.3606, 0.1394, 0.1394]
# Even split nucleotides
EVEN_DNA_PROBABILITIES = [0.25, 0.25, 0.25, 0.25]

SPECIES_NAME_DICT = {"hsa":"Homo sapiens", "pfa":"Plasmodium falciparum",
                     "sco":"Streptomyces coelicolor", "even":"50% GC content"}
RAND_SEQ_CLASS_LABEL = 0
RAND_IN_READ_STR_LABEL = 1
REPEAT_READ_LABEL = 2
FLANKING_READ_LABEL = 3

def compress(string):
    data = bytes(string, 'utf-8')
    # 1=Z_BEST_SPEED, -1=Z_DEFAULT_COMPRESSION, 9=Z_BEST_COMPRESSION
    s_out = zlib.compress(data, level=1)
    compression = (len(s_out)/len(data))
    return compression

def get_sequence_by_coord_local(selected_chromo, selected_start_pos,
                                selected_end_pos, selected_strand, reader):
    # TODO: strand selection
    return reader[selected_chromo][selected_start_pos:selected_end_pos].seq

def requests_retry_session(
    retries=3,
    backoff_factor=0.3,
    status_forcelist=(500, 502, 504),
    session=None,):
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

def get_sequence_by_coord(selected_chromo, selected_start_pos, selected_end_pos, selected_strand):
    server = "http://rest.ensembl.org"
    ext = "/sequence/region/human/"+selected_chromo+":"+str(selected_start_pos)+".."+str(selected_end_pos)+":"+ ("1" if selected_strand == "+" else "-1")+ "?coord_system_version=GRCh37;mask=soft"
    s = requests.Session()
    r = requests_retry_session(session=s).get(server + ext, headers={"Content-Type": "text/plain"})
    if not r.ok:
        print(r)
        r.raise_for_status()
        print("ERROR! Getter failed.")
        sys.exit(-2)
    else:
        return r.text.upper()

def get_random_DNA_string(genome_dict, read_len=150, reader=None):
    ''' Sample a random DNA sequence from hg19.

    :param read_len: the length of the DNA sample to return
    :return:
    '''
    sample_chr = random.choice(genome_dict)
    start_location = random.randint(1, int(sample_chr['length'])-read_len-1)
    strand = random.choice(["+", "-"])
    location_string = sample_chr['name']+":"+str(start_location)+"-"+ \
                      str(start_location+read_len)+":"+strand
    if reader:
        sample_seq = get_sequence_by_coord_local(sample_chr['name'], start_location,
                                                 (start_location+read_len), strand,
                                                 reader)
        if len(sample_seq) != read_len:
            print(len(sample_seq))
            print(location_string)
        assert len(sample_seq) == read_len
    else:
        sample_seq = get_sequence_by_coord(sample_chr['name'], start_location,
                                           (start_location+read_len-1), strand)
    return location_string, sample_seq


def get_random_pseudoDNA_string(probabilities, read_len=150):
    ''' Generate a random DNA-like string that has the same base distribution
    as human DNA (ie, CG proportion of 40.87% per Piovesan et al (2019), BMC
    Res. Notes).
    :param read_len: the length of the DNA-like string to return
    :return:
    '''
    random_seq = random.choices(DNA_ALPHABET, weights=probabilities,
                                k=read_len)
    return "".join(random_seq)


def get_repeat_string(repeat_length, repeat_motif=None):
    ''' Get a repeat length '''
    repeat_string = repeat_motif * int(repeat_length / len(repeat_motif))
    repeat_string += repeat_motif[0:repeat_length % len(repeat_motif)]
    return repeat_string


def get_DNA_random_samples(read_len=150, samplelimit=100000):
    if REFERENCE_FILE_PATH != "ENSEMBL":
        genome_dict_path = Path(REFERENCE_FILE_PATH).with_suffix(".superSTR.dict")
        if not genome_dict_path.exists():
            print("Generating contig length dictionary. This may take some time.")
            len_dict = defaultdict(int)
            with open(REFERENCE_FILE_PATH, 'rt') as in_ref:
                current_header = in_ref.readline().strip().replace(">", "")
                for line in in_ref:
                    if ">" in line:
                        current_header = line.strip().replace(">", "")
                    else:
                        len_dict[current_header] += len(line.strip())
            with open(genome_dict_path, 'wt') as genome_file:
                for key in len_dict:
                    genome_file.write(str(key)+"\t"+str(len_dict[key])+"\n")
        with open(genome_dict_path, 'rt') as genome_file:
            genome_dict = list(csv.DictReader(genome_file, fieldnames=["name", "length"], delimiter='\t'))
        reader = Fasta(REFERENCE_FILE_PATH)
    else:
        print("You have selected ENSEMBL download of human sequence. Depending on your internet connection this can take significant time; you may want to download this sequence and use the --reference_path option instead.")
        if not Path("hg19.chrom.sizes").exists():
            print("Downloading UCSC contig length dictionary...")
            url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes"
            response = requests.get(url)
            open("hg19.chrom.sizes", "wb").write(response.content)
        with open("hg19.chrom.sizes", 'rt') as genome_file:
            genome_dict = csv.DictReader(genome_file, fieldnames=["name", "length"], delimiter='\t')
            sublist = []
            for entry in genome_dict:
                key = entry["name"]
                if "chrUn" in key or "gl" in key or "hap" in key or "ctg" in key: continue
                if key == "chrM":
                    entry["name"] = "chrMT"
                sublist.append(entry)
            genome_dict = sublist
        reader = None
    compression_ratios = []
    sample_count = 0
    sublist = []
    for g_dict in genome_dict:
        if int(g_dict["length"]) <= read_len:
            print(g_dict["name"], "< read length of ", str(read_len), "dropping")
        else:
            sublist.append(g_dict)
    genome_dict = sublist
    with open("genome_data.txt", 'wt', buffering=1) as genome_dump:
        print("Obtaining sequence...")
        while sample_count < samplelimit:
            location_string, sample_seq = get_random_DNA_string(genome_dict,
                                                                read_len=read_len,
                                                                reader=reader)
            if sample_seq.count("N") <= 10:
                compression_ratios.append(compress(sample_seq))
                genome_dump.write(location_string)
                genome_dump.write("\t")
                genome_dump.write(sample_seq)
                genome_dump.write("\n")
                sample_count += 1
    return compression_ratios

def get_repeat_read(distribution, read_len=150, repeat_motif=None, max_repeat_period=15):
    '''Generate a read that contains only repetitive sequence.

    :param read_len (int): The read length; 150nt default.
    :param repeat_motif (str): Allows manual specification of repeat motif.
    :param max_repeat_period (int): Limit length of generated repeat.

    :return:
    '''
    # Get the period (ie, the k for the repeat k-mer) length.
    repeat_period = random.randint(1, max_repeat_period)
    # Generate a random repeat motif.
    if not repeat_motif:
        repeat_motif = "".join(random.choices(DNA_ALPHABET, weights=distribution, k=repeat_period))
    read = get_repeat_string(read_len, repeat_motif)
    descriptor_string = repeat_motif + "_" + str(read_len) \
                        + "bp_in_repeat_read"
    return descriptor_string, read


def get_random_in_read_str(distribution, read_len=150, repeat_motif="ACG",
                           max_repeat_period=None, min_number_repeats=25):
    '''Generate an in-read repeat. This is a short repetitive sequence flanked
    by random sequence on both sides within the read; ie, RAND-REPEAT-RAND.
    Such repeats are relatively short.

    :param read_len (int): The read length; 150nt default.
    :param repeat_motif (str): Allows manual specification of repeat motif.
    :param max_repeat_period (int): Limit length of generated repeat.
    :param min_number_repeats (int) : The minimum number of repeats of a motif.
    :param min_rand_seq_length (int): The minimum length of flanking random
        sequence

    :return:
    '''
    # The maximum amount of the in-read repeat that can be repetitive is 80%
    maximum_repeat_size = math.floor(read_len*0.8)
    if not max_repeat_period:
        max_repeat_period = math.floor(maximum_repeat_size/min_number_repeats)
    #print(str(maximum_repeat_size)+"/"+str(min_number_repeats) + " = " + str(max_repeat_period))
    # Get the period (ie, the k for the repeat k-mer) length.
    repeat_period = random.randint(2, max_repeat_period)
    # Generate a random repeat motif.
    if not repeat_motif:
        repeat_motif = "".join(random.choices(DNA_ALPHABET,
                                              weights=distribution,
                                              k=repeat_period))
    # Generate minimum and maximum repeat substring sizes
    minimum_repeat_size = repeat_period*min_number_repeats
    # Repeat length is random within the [minimum_repeat_size,
    # maximum_repeat_size] interval
    repeat_length = random.randint(minimum_repeat_size, maximum_repeat_size)
    avail_length = read_len-repeat_length
    flank_1_length = random.randint(3, avail_length)
    flank_2_length = avail_length-flank_1_length
    # # Pick a start position for the repeat.
    # startpos = random.randint(0, read_len - repeat_length)
    # Generate the lead substring; random sequence
    read = get_random_pseudoDNA_string(distribution, read_len=flank_1_length)
    # Generate the repeat substring; repeated STRs
    read += get_repeat_string(repeat_length, repeat_motif)
    # Generate the trailing substring; random sequence
    read += get_random_pseudoDNA_string(distribution, read_len=flank_2_length)
    # Generate descriptor string.
    descriptor_string = repeat_motif + "_" + str(repeat_length) \
                        + "repeat_in_read"
    return descriptor_string, read


def get_flanking_read_str(distribution, read_len=150, repeat_motif="ACG", left_assign=None,
                          max_repeat_period=None, min_number_repeats=10):
    '''Generate a flanking read. This is a short repetitive sequence joined to
    a random sequence; either RAND-REPEAT or REPEAT-RAND.

    :param read_len (int): The read length; 150nt default.
    :param repeat_motif (str): Allows manual specification of repeat motif.
    :param left_assign (bool): Allows manual specification of where the repeat
        is joined to the random sequence.
    :param max_repeat_period (int): Limit length of generated repeat.

    :return:
    '''
    maximum_repeat_size = math.floor(read_len * 0.8)
    # Get the period (ie, the k for the repeat k-mer) length.
    if not max_repeat_period:
        max_repeat_period = math.floor(maximum_repeat_size/min_number_repeats)
    repeat_period = random.randint(2, max_repeat_period)
    minimum_repeat_size = repeat_period*min_number_repeats
    # Get length of repetitive sequence.
    repeat_length = random.randint(minimum_repeat_size, maximum_repeat_size)
    nonrepeat_length = read_len - repeat_length
    # Generate a random repeat motif.
    if not repeat_motif:
        repeat_motif = "".join(random.choices(DNA_ALPHABET,
                                              weights=distribution,
                                              k=repeat_period))
    # Get flanking sequence string.
    seq_str = get_random_pseudoDNA_string(distribution, nonrepeat_length)
    # Get repeat sequence string.
    repeat_str = get_repeat_string(repeat_length, repeat_motif)
    # If we're not manually specifying the ordering, we randomly choose.
    if left_assign is None:
        left_assign = random.choice((True, False))
    if left_assign:
        # Read is sequence-repeat
        read = seq_str + repeat_str
        descriptor_string = repeat_motif + "_" + str(repeat_length) \
                            + "_seq_repeat_read"
    else:
        # Read is repeat-sequence
        read = repeat_str + seq_str
        descriptor_string = repeat_motif + "_" + str(repeat_length) \
                            + "_repeat_seq_read"
    return descriptor_string, read

def generate_pandas_sampleset(distribution, read_len=150):
    ratio_list = []
    label_list = []
    retain_list = []
    for i in range(0, 50000):
        if read_len > 25*3:
            ratio_list.append(
                compress(get_random_in_read_str(distribution, read_len)[1]))
            label_list.append(RAND_IN_READ_STR_LABEL)
            retain_list.append(1)
        if read_len > 10*3:
            ratio_list.append(
                compress(get_flanking_read_str(distribution, read_len)[1]))
            label_list.append(FLANKING_READ_LABEL)
            retain_list.append(1)
        ratio_list.append(compress(get_random_pseudoDNA_string(distribution, read_len)))
        label_list.append(RAND_SEQ_CLASS_LABEL)
        retain_list.append(0)
        ratio_list.append(compress(get_repeat_read(distribution, read_len)[1]))
        label_list.append(REPEAT_READ_LABEL)
        retain_list.append(1)
    for i in range(0, 50000*4):
        ratio_list.append(compress(get_random_pseudoDNA_string(distribution, read_len)))
        label_list.append(RAND_SEQ_CLASS_LABEL)
        retain_list.append(0)
    d = {"ratio":ratio_list, "class":label_list, "retain":retain_list}
    return pandas.DataFrame(d)

def render_distributions(random_sequence_ratios, repeat_in_read_ratios,
                         read_in_repeat_ratios, flanking_read_ratios, species,
                         dna_samples=None, read_len=150, threshold=0.55,
                         render_output_path="plots/"):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.distplot(repeat_in_read_ratios, hist=False, rug=False,
                 label="containing read")
    sns.distplot(read_in_repeat_ratios, hist=False, rug=False,
                 label="repeat read")
    sns.distplot(flanking_read_ratios, hist=False, rug=False,
                 label="flanking read")
    sns.distplot(random_sequence_ratios, hist=False, rug=False,
                 label="pseudorandom DNA/random read", color="black")
    if dna_samples:
        sns.distplot(dna_samples, hist=False, rug=False, label="DNA samples",
                     color="red")
    plt.title(""+str(read_len)+"nt read length")

    plt.xlabel("Compression ratio (compressed/uncompressed)")
    plt.ylabel("Density")
    plt.axvline(x=threshold, linestyle='--', color="grey")
    plt.savefig(render_output_path.joinpath(species+"_compressibility_classes_"+str(read_len)+"nt.png"), bbox_inches="tight")
    plt.close()

def get_random_quality_string(read_len):
    string_list = []
    for i in range(0, read_len):
        string_list.append(str(chr(random.randint(34, 125))))
    return "".join(string_list)

def generate_sampleset(distribution, read_len=150):
    if read_len > 30*3:
        repeat_in_read_ratios = [compress(
            get_random_in_read_str(distribution, read_len=read_len)[1])
                                 for i in range(0, 100000)]
    random_sequence_ratios = [compress(get_random_pseudoDNA_string(distribution, read_len=read_len))
                              for i in range(0, 100000)]
    read_in_repeat_ratios = [compress(get_repeat_read(distribution, read_len=read_len)[1])
                             for i in range(0, 100000)]
    flanking_read_ratios = [compress(get_flanking_read_str(distribution, read_len=read_len)[1])
                            for i in range(0, 100000)]
    return random_sequence_ratios, repeat_in_read_ratios, \
           read_in_repeat_ratios, flanking_read_ratios

def write_choice(choice_str, read_idx, read_len, file_handler, distribution):
    if choice_str == "random":
        file_handler.write("@rand_" + str(read_idx) + "\n" + get_random_pseudoDNA_string(distribution, read_len=read_len))
        file_handler.write("\n+\n"+get_random_quality_string(read_len)+"\n")
    elif choice_str == "rir":
        file_handler.write("@randir_" + str(read_idx) + "\n" + get_random_in_read_str(distribution, read_len=read_len)[1])
        file_handler.write("\n+\n" + get_random_quality_string(read_len) + "\n")
    elif choice_str == "pure":
        file_handler.write("@rir_" + str(read_idx) + "\n" + get_repeat_read(distribution, read_len=read_len)[1])
        file_handler.write("\n+\n" + get_random_quality_string(read_len) + "\n")
    elif choice_str == "flanking":
        file_handler.write("@flank_" + str(read_idx) + "\n" + get_flanking_read_str(distribution, read_len=read_len)[1])
        file_handler.write("\n+\n" + get_random_quality_string(read_len) + "\n")
    else:
        print("ERROR")
        import sys
        sys.exit(-1)

def generate_and_save_simulated_samples(distribution, output_root, species):
    out_root_folder = output_root.joinpath("simulations/")
    out_root_folder.mkdir(parents=True, exist_ok=True)
    out_random_path = out_root_folder.joinpath(species+"_random_seq.fa.gz")
    random_sequence_ratios = []
    with gzip.open(out_random_path, 'wt') as random_seq_file:
        for i in range(0, sample_size):
            seq = get_random_pseudoDNA_string(distribution,
                                              read_len=master_read_len)
            random_seq_file.write(">rand_" + str(i) + "\n")
            random_seq_file.write(seq)
            random_seq_file.write("\n")
            random_sequence_ratios.append(compress(seq))
    repeat_in_read_ratios = []
    out_repeat_in_read_path = out_root_folder.joinpath(species+"_repeat_in_read.fa.gz")
    with gzip.open(out_repeat_in_read_path, 'wt') as repeat_seq_file:
        for i in range(0, sample_size):
            detail_string, seq = get_random_in_read_str(distribution,
                                                        read_len=master_read_len)
            repeat_seq_file.write(
                ">repeat_in_read_" + str(i) + "_" + detail_string + "\n")
            repeat_seq_file.write(seq)
            repeat_seq_file.write("\n")
            repeat_in_read_ratios.append(compress(seq))
    read_in_repeat_ratios = []
    out_read_in_repeat_path = out_root_folder.joinpath(species+"_read_in_repeat.fa.gz")
    with gzip.open(out_read_in_repeat_path, 'wt') as repeatread_seq_file:
        for i in range(0, sample_size):
            detail_string, seq = get_repeat_read(distribution,
                                                 read_len=master_read_len)
            repeatread_seq_file.write(
                ">read_in_repeat_" + str(i) + "_" + detail_string + "\n")
            repeatread_seq_file.write(seq)
            repeatread_seq_file.write("\n")
            read_in_repeat_ratios.append(compress(seq))
    flanking_read_ratios = []
    out_flank_path = out_root_folder.joinpath(species+"_flanking.fa.gz")
    with gzip.open(out_flank_path, 'wt') as flank_seq_file:
        for i in range(0, sample_size):
            detail_string, seq = get_flanking_read_str(distribution,
                                                       read_len=master_read_len)
            flank_seq_file.write(
                ">flanking_" + str(i) + "_" + detail_string + "\n")
            flank_seq_file.write(seq)
            flank_seq_file.write("\n")
            flanking_read_ratios.append(compress(seq))

if __name__ == "__main__":
    # Argparser for command line arguments.
    import argparse
    parser = argparse.ArgumentParser()
    mode_group = parser.add_argument_group(title="Mode selection (one argument required)")
    mode_group.add_argument("--simulation", dest="simulation",
                            action="store_true",
                            help="Generate (and analyse) simulated data.")
    mode_group.add_argument("--performance", action="store_true",
                            dest="performance", help="Generate PRC and ROC "
                                                     "curves for the "
                                                     "compressor.")
    required_group = parser.add_argument_group(title="Required arguments")
    required_group.add_argument("--output", action="store", required=True, dest="output_path",
                        help="Output path for files.")
    parser.add_argument("--sample", action="store_true", dest="sample",
                        help="Use sampled DNA sequence, not pseudorandom. Requires either --reference_path or --ensembl flags set; the latter requires an internet connection.")
    parser.add_argument("-r", action="store", dest="read_len",
                        help="Read length for simulation (default: 150nt)")
    parser.add_argument("-t", action="store", dest="threshold",
                        help="Compressibility threshold.")
    parser.add_argument("--reference_path", action="store", dest="reference_path",
                        help="Path to genome reference for sampling; must be uncompressed FASTA")
    parser.add_argument("--write_sims", action="store_true", dest="write_sims",
                        help="Write simulated reads to disk; requires output "
                             "path.")

    parser.add_argument("--sample_size", action="store", dest="sample_size",
                        help="Number of sequences to generate/sample.")
    parser.add_argument("--hsa", action="store_true", dest="hsa",
                        default=False,
                        help="Simulate with H. sapiens genomic GC content.")
    parser.add_argument("--sco", action="store_true", dest="sco",
                        default=False,
                        help="Simulate with S. coelicolor genomic GC content.")
    parser.add_argument("--pfa", action="store_true", dest="pfa",
                        default=False,
                        help="Simulate with P. falciparum genomic GC content.")
    parser.add_argument("--even", action="store_true", dest="even",
                        default=False,
                        help="Simulate with even probabilities for each "
                             "nucleotide")
    parser.add_argument("--all", action="store_true", dest="all", default=False,
                        help="Simulate for all four GC contents.")

    parser.add_argument("--ensembl", action="store_true", dest="ensembl", default=False,
                        help="Use randomly selected ENSMEBL hg19 sequence instead of pseudorandom genomic sequence for background.")
    if len(sys.argv) == 1:
        print("ERROR: No options specified.", file=sys.stderr)
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    if not (args.simulation or args.performance):
        print("ERROR: Select at least one of --simulation and --performance.", file=sys.stderr)
        parser.print_help(sys.stderr)
        sys.exit(1)
    # Parse variables from input, otherwise use sensible defaults.
    master_read_len = args.read_len if args.read_len else 150
    sample_size = args.sample_size if args.sample_size else 150000
    if args.all:
        species_rundict = {"hsa":True, "sco":True, "pfa":True, "even":True}
    else:
        species_rundict = {"hsa": args.hsa, "sco": args.sco, "pfa": args.pfa,
                           "even": args.even}
    dist_dict = {"hsa":HSA_DNA_PROBABILITIES, "sco":SCO_DNA_PROBABILITIES,
                 "pfa":PFA_DNA_PROBABILITIES, "even":EVEN_DNA_PROBABILITIES}
    if args.simulation:
        if not args.threshold:
            print("ERROR: --simulation mode requires -t to be set.", file=sys.stderr)
            parser.print_help(sys.stderr)
            sys.exit(1)
        threshold = float(args.threshold)
        for species in species_rundict:
            dist = dist_dict[species]
            if not species_rundict[species]:
                continue
            # Generate simulated files:
            if args.write_sims:
                random_sequence_ratios, repeat_in_read_ratios, \
                read_in_repeat_ratios, flanking_read_ratios = \
                    generate_and_save_simulated_samples(dist,
                                                        Path(args.output_path),
                                                        species)
            else:
                random_sequence_ratios, repeat_in_read_ratios, \
                read_in_repeat_ratios, flanking_read_ratios = \
                    generate_sampleset(dist, read_len=master_read_len)
            if args.sample:
                if args.reference_path:
                    REFERENCE_FILE_PATH = args.reference_path
                elif args.ensembl:
                    REFERENCE_FILE_PATH = "ENSEMBL"
                else:
                    print("ERROR: Sampling mode requires one of --reference_path or --ensembl to be set.", file=sys.stderr)
                    parser.print_help(sys.stderr)
                    sys.exit(1)
                random_dna_ratios = get_DNA_random_samples(read_len=master_read_len)
                render_distributions(random_dna_ratios, repeat_in_read_ratios,
                                     read_in_repeat_ratios, flanking_read_ratios,
                                     species, read_len=master_read_len,
                                     threshold=threshold,
                                     render_output_path=Path(args.output_path))
            else:
                render_distributions(random_sequence_ratios,
                                     repeat_in_read_ratios,
                                     read_in_repeat_ratios,
                                     flanking_read_ratios, species,
                                     read_len=master_read_len,
                                     threshold=threshold,
                                     render_output_path=Path(args.output_path))
    elif args.performance:
        for species in species_rundict:
            if not species_rundict[species]:
                continue
            print(SPECIES_NAME_DICT[species])
            print("-"*80)
            print("Read length\tROC\tPrecision+Recall")
            if args.read_len and str.isdigit(args.read_len):
                proc_list = [int(args.read_len)]
            else:
                proc_list = [16, 36, 75, 150, 200, 250, 300]
            for rlen in proc_list:
                data_frame = generate_pandas_sampleset(dist_dict[species], read_len=rlen)
                retain_nparray = data_frame['retain'].to_numpy()
                ratio_nparray = data_frame['ratio'].to_numpy()
                fpr, tpr, thresholds = metrics.roc_curve(retain_nparray, ratio_nparray,
                                                         pos_label=0)
                optimal = numpy.argmax(tpr-fpr)
                roc_threshold = thresholds[optimal]
                roc_auc = auc(fpr, tpr)
                import matplotlib.pyplot as plt
                plt.figure()
                lw = 2
                plt.plot(fpr, tpr, color='darkorange',
                         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
                plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC - 400,000 simulated reads\n50% retain; min ')
                plt.legend(loc="lower right")
                render_path = Path(args.output_path)
                plt.savefig(render_path.joinpath(args.output_path, species+"_ROC"+str(rlen)+"nt.png"), bbox_inches="tight")
                plt.close()
                precision, recall, thresholds = precision_recall_curve(retain_nparray,
                                                                       ratio_nparray,
                                                                       pos_label=0)
                optimal = numpy.argmax(precision + recall)
                prc_threshold = thresholds[optimal]
                print(str(rlen) + "\t" + str(roc_threshold) + "\t" + str(prc_threshold))
                plt.figure()
                lw = 2
                plt.plot(recall, precision, color='darkorange',
                         lw=lw, label='PRC curve (asrea = %0.2f)' % roc_auc)
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.ylabel('Precision')
                plt.xlabel('Recall')
                plt.title(
                    'PRC - gzip compression of 400,000 simulated reads\n(200,000 '
                    'retain, 200,000 reject)')
                plt.legend(loc="lower right")
                plt.savefig(Path(args.output_path).joinpath(args.output_path, species+"_PRC_"+str(rlen)+"nt.png"), bbox_inches="tight")
                plt.close()