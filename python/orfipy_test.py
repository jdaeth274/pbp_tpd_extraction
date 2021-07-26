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

toc_importing = time.perf_counter()

print("Importing packages took %s" % (toc_importing - tic_importing))


#############
# FUNCTIONS #
#############

def gene_name_search(attribute_col, gene_name):
    gene_finder = "gene=" + gene_name + ";"
    search_res = re.findall(gene_finder, attribute_col)
    if search_res == []:
        search_res = "No Gene Name in GFF"
    else:
        search_res = str(search_res)
        search_res = re.findall('=.*?;', search_res)
        search_res = str(search_res)

        search_res = search_res[3:-3]

    return search_res


def fasta_file_searcher(base_name, fasta_list):
    bassy_name = os.path.basename(base_name)

    bassy_name = re.split("\.", bassy_name, maxsplit=2)[0]

    bassy_name = bassy_name + ".contigs_velvet.fa"

    fasta_path = list(filter(lambda x: bassy_name in x, fasta_list))

    return fasta_path[0]


def get_options():
    purpose = ''' This is a script to input a set of fastas, extract genes of
    interest and then search these against a database for certain proteins.
    Usage: python pen_checker_cdc.py <input_gffs> <fasta_list> <gene_name> <backup_gene_name> <gene_alignment> <target_length> <output_name>'''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='pen_checker_cdc.py')

    parser.add_argument('--seqs', required=True, help='List of seqeuence files (GFF or FASTA) (required)', type=str)
    parser.add_argument('--gene', required=True, help='Specify PBP being analysed  (required)', type=str,
                        choices=['pbp1a', 'pbp2b', 'pbp2x', 'dhfR', 'folP'])
    parser.add_argument('--cores', required=False, help='Number of cores to use for orfipy',
                        default=None)
    parser.add_argument('--aln', required=False, help='Gene alignment  (required)', type=str)
    parser.add_argument('--tlength', required=True, help='Target length of TPD domain  (required)', type=int)
    parser.add_argument('--tolerance', help='Deviation from expected length allowed', type=int, default=10)
    parser.add_argument('--data_dir', help='Location of the data directory of package', type=str, default="./data")
    parser.add_argument('--output', required=True, help='Prefix of output files  (required)', type=str)

    args = parser.parse_args()

    return args


def search_for_gene(ref_in, name, gene_length, tol, correct_length, gene_rower):
    gene_rowz = gene_rower
    if not correct_length:
        correct_length = False
        gene_row = ref_in['attributes'].str.contains(name)
        gene_row_indy = gene_row.where(gene_row == True)
        gene_row_indy = gene_row_indy.index[gene_row_indy == True].tolist()
        gene_rowz = ref_gff_tsv.iloc[gene_row_indy]
        if gene_rowz.empty == False:
            gene_len = [abs(int(gene_rowz.iloc[0, 4]) - int(gene_rowz.iloc[0, 3]))]

            overhang = [gene_length - tol, gene_length + tol]

            if overhang[0] <= gene_len[0] <= overhang[1]:
                correct_length = True
            else:
                correct_length = False
                sys.stderr.write('Found gene' + name + ' but wrong length: ' + str(gene_len[0]) + ', expected: ' + str(
                    gene_length) + '\n')

    return correct_length, gene_rowz


def get_aln_pos_from_ref(hmm_aln, pos, offset):
    pos = pos - offset
    ref_pos = pos + 1
    upstream_length = 0
    while_counter = 0
    while upstream_length < pos:
        old_upstream = upstream_length
        upstream_frag = hmm_aln[0, 0:ref_pos].seq

        upstream_length = len(upstream_frag) - upstream_frag.count('.')
        ref_pos = ref_pos + 1
        while_counter += 1
        if while_counter > 10:
            if (upstream_length - old_upstream) == 0:
                return False

    return (ref_pos - 1)

class seqqer():
    def __init__(self):
        self.seq = None


def hmm_search_for_gene(fasta, gene, aa_dir_name, data_dir, gene_leng, tolerance, num_cores):
    # data structure
    readingframes = []

    # from https://www.biostars.org/p/183616/
    # load sequence

    tic_aa_creator = time.perf_counter()
    print("Creating aa")

    ## run orfipy ORF finder on the fasta file

    orfipy_run = "orfipy " + fasta + " --pep tmp_out_pep.fa --outdir tmp_orfi_out --min 50 --procs " + str(num_cores)

    subprocess.run(orfipy_run, shell=True)

    aa_file = "./tmp_orfi_out/tmp_out_pep.fa"
    #
    #
    # for record in SeqIO.parse(fasta, "fasta"):
    #
    #     # process file name expecting no extra .
    #     gff_base = os.path.basename(fasta)
    #     gff_base = re.sub("\..*[a-z,A-Z,0-9].*$", "", gff_base)
    #     aa_base_name = aa_dir_name + "/" + gff_base
    #
    #     # Create three reading frames in forward direction, offset 0, 1, 2
    #     record_rc = record.reverse_complement()
    #     readingframes.extend(
    #         [Seq.translate(record.seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in
    #          range(3)])
    #     readingframes.extend(
    #         [Seq.translate(record_rc.seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in
    #          range(3)])
    #
    #     results = []
    #     for frame in readingframes:
    #         for peptide in frame.split('*'):  # Split translation over stopcodons
    #             if len(peptide) > 30:
    #                 results.append(peptide)
    #
    #     # Write length and translation to file
    #     # Use PotentialORFs.txt as output, can be changed
    #     with open(aa_base_name + '.aa', 'w') as output:
    #         for n, peptide in enumerate(results):
    #             output.write(">{}\n{}\n".format('peptide_' + str(n), peptide))

    toc_aa_creator = time.perf_counter()
    print("Running HMM on aa")
    tic_hmm_run = time.perf_counter()
    # run HMM
    gff_base = os.path.basename(fasta)
    gff_base = re.sub("\..*[a-z,A-Z,0-9].*$", "", gff_base)
    aa_base_name = aa_dir_name + "/" + gff_base
    hmm_output_fn = aa_base_name
    hmm_proc_out = 1
    while hmm_proc_out != 0:
        hmm_proc_out = subprocess.check_call(
            'hmmsearch ' + data_dir + '/' + gene + '.hmm ' + aa_file + ' > ' + hmm_output_fn, shell=True)

    subprocess.run("rm -r tmp_orfi_out", shell=True)

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
    hmm_offset = hmm_best_hsp.query_start
    hmm_aln = hmm_best_hsp.aln
    subprocess.call('rm ' + hmm_output_fn, shell=True)

    toc_hmm_run = time.perf_counter()
    print("Printing out the resistance")
    seq = None
    if gene == 'folP':
        aln_start = get_aln_pos_from_ref(hmm_aln, 57, hmm_offset)
        aln_end = get_aln_pos_from_ref(hmm_aln, 67, hmm_offset)
        if isinstance(aln_start, bool) or isinstance(aln_end, bool):
            print(gff_base + '\tMatch not large enough to folP')
            print("Took this long for ORF finder: %s" % (toc_aa_creator - tic_aa_creator))
            print("Took this long for HMM run: %s" % (toc_hmm_run - tic_hmm_run))
            status = "NA"
            return status, gff_base
        hmm_match = hmm_aln[0, aln_start:aln_end].seq
        query_match = hmm_aln[1, aln_start:aln_end].seq
        gap_count = hmm_match.count('.') - query_match.count('.')

        if gap_count > 0:
            print(gff_base + '\tSulphamethoxazole resistant')
            status = "R"
        else:
            print(gff_base + '\tSulphamethoxazole sensitive')
            status = "S"

    elif gene == 'dhfR':
        print("Getting aln pos")
        aln_start = get_aln_pos_from_ref(hmm_aln, 98, hmm_offset)
        if isinstance(aln_start, bool):
            print(gff_base + '\tMatch not large enough to dhfR')
            print("Took this long for ORF finder: %s" % (toc_aa_creator - tic_aa_creator))
            print("Took this long for HMM run: %s" % (toc_hmm_run - tic_hmm_run))
            status = "NA"
            return status, gff_base
        hmm_match = hmm_aln[0, aln_start:(aln_start + 1)].seq
        query_match = hmm_aln[1, aln_start:(aln_start + 1)].seq
        print("Printing out results")
        if query_match.upper() == "L":
            print(gff_base + '\tTrimethoprim resistant (' + query_match.upper() + ')')
            status = "R"
        elif query_match.upper() == "I":
            print(gff_base + '\tTrimethoprim sensitive (' + query_match.upper() + ')')
            status = "S"
        else:
            print(gff_base + '\tTrimethoprim unknown: ' + str(query_match).upper())
            status = "NA"

    elif gene in ['pbp1a','pbp2b','pbp2x']:
        aa_seq = hmm_best_hsp.hit
        if len(aa_seq) < gene_leng - tolerance:
            seq = seqqer()
            print("Gene length not long enough has %s needs %s" % (len(aa_seq), gene_leng))
        else:
            seq = aa_seq
        status = None

    print("Took this long for ORF finder: %s" % (toc_aa_creator - tic_aa_creator))
    print("Took this long for HMM run: %s" % (toc_hmm_run - tic_hmm_run))

    return status, gff_base, seq


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


def process_pbp(isolate, protein_fasta, tpd_start, tpd_end, tpd_lab, results_csv, k):
    protein_rec = list(SeqIO.parse(protein_fasta, format="fasta"))
    prot = protein_rec[0].seq

    tpd_string = str(prot[tpd_start: tpd_end])

    tpd_string = Seq(tpd_string)  # , generic_protein)

    prot_tpd_id = isolate + "_" + tpd_lab + "_TPD"
    tpd_protein_string = SeqIO.SeqRecord(tpd_string, id=prot_tpd_id)

    prot_file = isolate + "_" + tpd_lab + ".prot"
    with open(prot_file, "w+") as output_handle:
        SeqIO.write(tpd_protein_string, output_handle, "fasta")

    top_id = str(results_csv.iloc[0, 1])
    top_id = top_id[2:]

    if bassio_nameo != "22841_3#15.contigs_velvet.fa.gff":

        rm_command = "rm " + protein_fasta + " " + protein_csv
        os.system(rm_command)
    else:
        print(sstart, send)
        print(tpd_start, tpd_end)

    print("Generating CSV: %s%%" % round((k / len(gff_lines) * 100)))


def contig_number_getter(contig_name):
    if "|" in contig_name:
        contig_num = contig_name[-3:]
        contig_num = int(contig_num)
    elif "NODE" in contig_name:
        contig_num = re.sub("_length.*$", "", contig_name)
        contig_num = re.sub("NODE_", "", contig_num)
        contig_num = int(contig_num)
    else:
        contig_num = re.sub("^.*[0-9]\.", "", contig_name)
        contig_num = int(contig_num)

    return contig_num


#########
# BEGIN #
#########

if __name__ == '__main__':
    tic_setup = time.perf_counter()
    print("Beginning setup")
    pandas.set_option('display.max_columns', 500)

    # parse command line
    files_for_input = get_options()

    ###############################################################################
    ## Lets go through the gff list in a for loop, we'll find the position of the #
    ## gene of interest and then from there extract this sequence and compare it ##
    ## to the gene alignment to test which is the right gene ######################
    ###############################################################################
    gff_files = open(files_for_input.seqs, "r")
    gff_lines = gff_files.read().splitlines()

    if gff_lines[0].endswith(".gff"):
        fasta_lines = []
        for gff in gff_lines:
            fasta_lines.append(extract_fasta_from_gff(gff))
    else:
        fasta_lines = gff_lines

    if files_for_input.cores is None:
        num_cores = mp.cpu_count() - 1
    else:
        num_cores = files_for_input.cores

    gene_name = files_for_input.gene
    # input_gene_names = files_for_input.gene.split(',')  # allow a list of alternative names
    # input_gene_names.append(gene_name)
    # input_gene_names.append(gene_name.lower())
    # all_gene_names = set(input_gene_names)

    gene_length = files_for_input.tlength


    aa_dir_name = "./" + re.sub(".csv", "", files_for_input.output) + "_aa_dir"
    if not os.path.exists(aa_dir_name):
        os.mkdir(aa_dir_name)
    all_gene_names = [files_for_input.gene]

    ###############################################################################
    ## So we've loaded up our collection of fastas and our gffs collection, at ####
    ## the moment these is lined up to work with contig files, as the gff annots ##
    ## have the start and end as the contig positions, this will need to be #######
    ## altered if we want to use just whole dna files.                          ###
    ## For now though we'll loop through the gff files, find the gene's contig ####
    ## and position, then using these we'll get the seq from the fasta list ,then #
    ## we'll convert these to protein and then blast against the reference align ##
    ## and from there pick the top blast hit ######################################
    ###############################################################################

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


        ref_gff_tsv = pandas.read_csv(gff_file, sep='\t',
                                      names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                             'attributes'],
                                      header=None)

        if ref_gff_tsv['attributes'].isnull().all():
            missing_isolates.append(bassio_nameo)
            missing_gff_isolates.append(bassio_nameo)
            if files_for_input.gene in ['pbp1a', 'pbp2b', 'pbp2x']:
                print(("This isn't a recognised GFF: " + bassio_nameo))
                continue

        gene_rower = None
        correct_length = False

        current_res, current_isolate,seq = hmm_search_for_gene(fasta_file, files_for_input.gene, aa_dir_name,
                                                           files_for_input.data_dir, gene_length, files_for_input.tolerance,
                                                               num_cores)

        if gene_name not in ['pbp1a', 'pbp1A', 'pbp2b', 'pbp2B', 'pbp2x', 'pbp2X']:
            gene_id.append(current_res)

            isolate_name.append(current_isolate)
            skip = True
        toc_iso_run = time.perf_counter()
        print("Took this long for isolate %s, overall: %s" % (current_isolate, (toc_iso_run - tic_iso_run)))
        if gene_name not in ['pbp1a','pbp1A','pbp2b','pbp2B','pbp2x','pbp2X']:
            continue


            # write protein sequence out for BLAST
        if seq.seq == None:
            continue
        isolate = re.split("\.", bassio_nameo)[0]
        with open((bassio_nameo + files_for_input.gene + ".fasta"), "w+") as output_handle:
            SeqIO.write(seq, output_handle, "fasta")

            #######################################################################
            ## Now we've got our protein string we'll use blast to search for the #
            ## closest hit in the cdc alignment file ##############################
            #######################################################################

        protein_fasta = bassio_nameo + files_for_input.gene + ".fasta"
        protein_csv = bassio_nameo + files_for_input.gene + ".csv"

        align_path = files_for_input.aln
        basename_aligno = os.path.basename(files_for_input.aln)
        basename_blastdb = basename_aligno + ".phr"

        if not os.path.isfile(basename_blastdb):
            db_command = "makeblastdb -in " + align_path + " -dbtype prot -out " + basename_aligno
            subprocess.call(db_command, shell=True)

        blasty_cmd = "blastp -query " + protein_fasta + " -db " + basename_aligno + " -out " + protein_csv + " -outfmt 10"
        subprocess.call(blasty_cmd, shell=True)

        #######################################################################
        ## Now we've got the blast results, lets extract the top hit and then #
        ## assign this isolates name as that particular allele type ###########
        #######################################################################

        results_csv = pandas.read_csv(protein_csv,
                                      names=['qid', 'sid', 'pid', 'align',
                                             'gap', 'mismatch', 'qstart', 'qend',
                                             'sstart', 'send', 'eval', 'bitscore'])

        results_csv = results_csv.sort_values(['bitscore'], ascending=False)
        if results_csv.empty:

            print("No Blast results for this isolate:", bassio_nameo)
            sys.exit(1)
            missing_isolates.append(bassio_nameo)
        else:
            top_res = results_csv.iloc[0]
            sstart = top_res.iloc[6] - 1
            send = top_res.iloc[7]

            if files_for_input.gene == "pbp1a":
                tpd_start = sstart
                tpd_end = send  # sstart + 252
                tpd_lab = "pbp1a"
            elif files_for_input.gene == "pbp2b":
                tpd_start = sstart
                tpd_end = send  # sstart + 277
                tpd_lab = "pbp2b"
            elif files_for_input.gene == "pbp2x":
                tpd_start = sstart  # sstart + 60
                tpd_end = send  # tpd_start + 299
                tpd_lab = "pbp2x"
                if isolate == "11511_7#11":
                    print(sstart, tpd_start, "Now for the ends \n", send, tpd_end)
                    print(top_res)
            else:
                sys.exit(
                    "For TPD extraction please use one of the three following gene names: pbp1a, pbp2b or pbp2x")

            isolate_name.append(isolate)
            top_pbp_hit = process_pbp(isolate, protein_fasta, tpd_start, tpd_end, tpd_lab, results_csv, k)
            gene_id.append(top_pbp_hit)

    out_df = pandas.DataFrame()
    out_df['isolate_id'] = pandas.Series(isolate_name)
    if files_for_input.gene in ['pbp1a', 'pbp2b', 'pbp2x']:
        out_df['top_id'] = pandas.Series(gene_id, index=out_df.index)
    else:
        out_df['Resistance'] = pandas.Series(gene_id, index=out_df.index)
        subprocess.call("rm -r *_aa_dir", shell=True)

    # rm_command = "rm " + protein_fasta + " " + protein_csv
    # os.system(rm_command)

    if len(missing_isolates) > 0:
        with open(("./missing_" + files_for_input.gene + "_ids.txt"), mode='wt', encoding='utf-8') as myfile:
            myfile.write('\n'.join(missing_isolates) + '\n')
    if len(missing_gff_isolates) > 0:
        with open(("./missing_gff_" + files_for_input.gene + "_ids.txt"), mode='wt', encoding='utf-8') as myfile:
            myfile.write('\n'.join(missing_gff_isolates) + '\n')
    out_df.to_csv(files_for_input.output,
                  index=False)













