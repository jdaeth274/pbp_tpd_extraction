print("Importing packages")
import time

tic_importing = time.perf_counter()
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SearchIO
import pandas
import re
import os
import sys
import argparse
import subprocess
import multiprocessing as mp

toc_importing = time.perf_counter()

print("Importing packages took %s" % (toc_importing - tic_importing))


#############
# FUNCTIONS #
#############

def get_options():
    purpose = ''' This is a script to extract the full lengths of the pbp sequences from 
    gff files. input: <gff_list> <data_dir> '''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='pen_checker_cdc.py')

    parser.add_argument('--gff', required=True, help='List of GFF files (required)', type=str)
    parser.add_argument('--data_dir', help='Location of the data directory of package', type=str, default="./data")
    parser.add_argument('--output', required=True, help='Prefix of output files  (required)', type=str)
    parser.add_argument('--cores', help='Number of cores available to run (default num_cores -1)', type=int, default=None)

    args = parser.parse_args()

    return args



class seqqer():
    def __init__(self):
        self.seq = None

def hmm_search_for_gene(fasta, data_dir, num_cores):
    # data structure
    readingframes = []

    # from https://www.biostars.org/p/183616/
    # load sequence

    tic_aa_creator = time.perf_counter()
    print("Creating aa")

    ## run orfipy ORF finder on the fasta file

    orfipy_run = "orfipy " + fasta + " --pep tmp_out_pep.fa --outdir tmp_orfi_out --min 50 --procs " + num_cores

    subprocess.run(orfipy_run, shell=True)

    toc_aa_creator = time.perf_counter()
    # print("Running HMM on aa")
    # tic_hmm_run = time.perf_counter()
    # # run HMM
    # gff_base = os.path.basename(fasta)
    # gff_base = re.sub("\..*[a-z,A-Z,0-9].*$", "", gff_base)
    # aa_base_name = aa_dir_name + "/" + gff_base
    # hmm_output_fn = aa_base_name
    # hmm_proc_out = 1
    # while hmm_proc_out != 0:
    #     hmm_proc_out = subprocess.check_call(
    #         'hmmsearch ' + data_dir + '/' + gene + '.hmm ' + aa_file + ' > ' + hmm_output_fn, shell=True)
    #
    # subprocess.run("rm -r tmp_orfi_out", shell=True)
    #
    # # parse HMM
    # print("Parsing HMM")
    # hmm_output = list(SearchIO.parse(hmm_output_fn, 'hmmer3-text'))
    # hmm_best_hsp_index = None
    # hmm_best_hsp_bitscore = 0
    # for i, hsp in enumerate(hmm_output[0]):
    #     if hsp.bitscore > hmm_best_hsp_bitscore:
    #         hmm_best_hsp_bitscore = hsp.bitscore
    #         hmm_best_hsp_index = i
    #
    # if hmm_best_hsp_index is None:
    #     return "Missing_HMM", gff_base
    #
    # hmm_best_hsp = hmm_output[0][hmm_best_hsp_index][0]
    # hmm_offset = hmm_best_hsp.query_start
    # hmm_aln = hmm_best_hsp.aln
    # subprocess.call('rm ' + hmm_output_fn, shell=True)
    #
    # toc_hmm_run = time.perf_counter()
    #
    #
    # aa_seq = hmm_best_hsp.hit
    # if len(aa_seq) < gene_leng - tolerance:
    #     seq = seqqer()
    #     print("Gene length not long enough has %s needs %s" % (len(aa_seq), gene_leng))
    # else:
    #     seq = aa_seq
    # status = None
    seq1a = hmm_run(fasta, data_dir, "pbp1a",720, 100)
    seq2b = hmm_run(fasta, data_dir, "pbp2b",681, 100)
    seq2x = hmm_run(fasta, data_dir, "pbp2x", 751, 100)
    subprocess.run("rm -r tmp_orfi_out", shell=True)

    print("Took this long for ORF finder: %s" % (toc_aa_creator - tic_aa_creator))


    return seq1a, seq2b, seq2x

def hmm_run(fasta, data_dir, gene,gene_leng, tolerance):
    print("Running HMM on aa %s" % gene)
    aa_file = "./tmp_orfi_out/tmp_out_pep.fa"
    tic_hmm_run = time.perf_counter()
    # run HMM
    gff_base = os.path.basename(fasta)
    gff_base = re.sub("\..*[a-z,A-Z,0-9].*$", "", gff_base)
    aa_base_name = gff_base
    hmm_output_fn = aa_base_name + "_hmmsearch_res.txt"
    hmm_proc_out = 1
    while hmm_proc_out != 0:
        hmm_proc_out = subprocess.check_call(
            'hmmsearch ' + data_dir  + gene + '.hmm ' + aa_file + ' > ' + hmm_output_fn, shell=True)



    # parse HMM
    print("Parsing HMM")
    hmm_output = list(SearchIO.parse(hmm_output_fn, 'hmmer3-text'))
    hmm_best_hsp_index = None
    hmm_best_hsp_bitscore = 0
    for i, hsp in enumerate(hmm_output[0]):
        if hsp.bitscore > hmm_best_hsp_bitscore:
            hmm_best_hsp_bitscore = hsp.bitscore
            hmm_best_hsp_index = i

    if hmm_best_hsp_index is None:
        return "Missing_HMM", gff_base

    hmm_best_hsp = hmm_output[0][hmm_best_hsp_index][0]
    subprocess.call('rm ' + hmm_output_fn, shell=True)

    toc_hmm_run = time.perf_counter()

    aa_seq = hmm_best_hsp.hit
    if len(aa_seq) < gene_leng - tolerance:
        seq = seqqer()
        print("Gene length not long enough has %s needs %s" % (len(aa_seq), gene_leng))
    else:
        seq = aa_seq
    status = None
    print("Took this long for HMM run: %s" % (toc_hmm_run - tic_hmm_run))

    return seq

def extract_fasta_from_gff(gff_fn):
    # check file name
    fa_fn = gff_fn
    if gff_fn.endswith('.gff'):
        fa_fn = fa_fn.replace('.gff', '.fa')
    else:
        sys.stderr.write('Does not appear to be a GFF file: ' + gff_fn + '\n')
        exit(1)

    # write output
    in_fa_region = False
    with open(gff_fn, 'r') as gff_file, open(fa_fn, 'w') as fa_file:
        for line in gff_file.readlines():
            if line.startswith('>'):
                in_fa_region = True
            if in_fa_region:
                fa_file.write(line)

    # return name
    return fa_fn




#########
# BEGIN #
#########

if __name__ == '__main__':
    tic_setup = time.perf_counter()
    print("Beginning setup")
    pandas.set_option('display.max_columns', 500)

    # parse command line
    files_for_input = get_options()

    ## Get number of cores available
    if files_for_input.cores is None:
        num_cores = mp.cpu_count() - 1
    else:
        num_cores = files_for_input.cores


    gff_files = open(files_for_input.gff, "r")
    gff_lines = gff_files.read().splitlines()

    fasta_lines = []
    for gff in gff_lines:
        fasta_lines.append(extract_fasta_from_gff(gff))
    # else:
    #     fasta_lines = open(files_for_input.fasta, "r")
    #     fasta_lines = fasta_lines.read().splitlines()

    #gene_name = files_for_input.gene
    # input_gene_names = files_for_input.gene.split(',')  # allow a list of alternative names
    # input_gene_names.append(gene_name)
    # input_gene_names.append(gene_name.lower())
    # all_gene_names = set(input_gene_names)

    #gene_length = files_for_input.tlength


    aa_dir_name = "./" + re.sub(".csv", "", files_for_input.output) + "_aa_dir"

    isolate_name = []
    gene_id = []
    missing_isolates = []
    missing_gff_isolates = []
    skip = False
    toc_setup = time.perf_counter()
    print("Took this long for initial set up: %s" % (toc_setup - tic_setup))
    print("Beginning iso run")
    for k, (gff_file, fasta_file) in enumerate(zip(gff_lines, fasta_lines)):
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("On isolate %s of %s" % (k, len(gff_lines)))
        tic_iso_run = time.perf_counter()
        bassio_nameo = os.path.basename(gff_file)
        isolate = re.split("\.", bassio_nameo)[0]

        gene_rower = None
        correct_length = False

        seq1a, seq2b, seq2x = hmm_search_for_gene(fasta_file, files_for_input.data_dir, num_cores)
        if seq1a.seq == None:
            print("pbp1a hit not long enough for this isolate")
        else:
            isolate = re.split("\.", bassio_nameo)[0]
            with open((isolate + "_pbp1a" + ".prot"), "w+") as output_handle:
                SeqIO.write(seq1a, output_handle, "fasta")
        if seq2b.seq == None:
            print("pbp2b hit not long enough for this isolate")
        else:
            isolate = re.split("\.", bassio_nameo)[0]
            with open((isolate + "_pbp2b" + ".prot"), "w+") as output_handle:
                SeqIO.write(seq2b, output_handle, "fasta")

        if seq2x.seq == None:
            print("pbp2x hit not long enough for this isolate")
        else:
            isolate = re.split("\.", bassio_nameo)[0]
            with open((isolate + "_pbp2x" + ".prot"), "w+") as output_handle:
                SeqIO.write(seq2x, output_handle, "fasta")

        pbp1a_cat = "ls -d $PWD/*pbp1a.prot > pbp1a_prot_list"
        pbp2b_cat = "ls -d $PWD/*pbp2b.prot > pbp2b_prot_list"
        pbp2x_cat = "ls -d $PWD/*pbp2x.prot > pbp2x_prot_list"

        subprocess.run(pbp1a_cat, shell=True)
        subprocess.run(pbp2b_cat, shell=True)
        subprocess.run(pbp2x_cat, shell=True)

        aa_df_creator = "Rscript --vanilla " + files_for_input.data_dir + "../R/full_seq_AA.R" + " pbp1a_prot_list pbp2b_prot_list pbp2x_prot_list " + files_for_input.output + ".csv"
        subprocess.run(aa_df_creator, shell=True)



