from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SearchIO
import pandas
import re
import os
import sys
import argparse
import subprocess
import time
#############
# FUNCTIONS #
#############

def gene_name_search(attribute_col, gene_name ):
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

    fasta_path = list(filter(lambda x:bassy_name in x,fasta_list ))

    return fasta_path[0]

def get_options():

    purpose = ''' This is a script to input a set of fastas, extract genes of
    interest and then search these against a database for certain proteins.
    Usage: python pen_checker_cdc.py <input_gffs> <fasta_list> <gene_name> <backup_gene_name> <gene_alignment> <target_length> <output_name>'''
    
    parser = argparse.ArgumentParser(description = purpose,
                                     prog='pen_checker_cdc.py')

    parser.add_argument('--gff', required=True, help='List of GFF files (required)', type = str)
    parser.add_argument('--fasta', required=False, help='List of FASTA files (ordered as GFF files)', type = str, default = None)
    parser.add_argument('--pbp', required=True, help='Specify PBP being analysed  (required)', type = str, choices=['pbp1a', 'pbp2b', 'pbp2x','dhfR','folP'])
    parser.add_argument('--gene', required=True, help='Comma-separated list of gene synonyms (required)', type = str)
    parser.add_argument('--aln', required=True, help='Gene alignment  (required)', type = str)
    parser.add_argument('--tlength', required=True, help='Target length of TPD domain  (required)', type = int)
    parser.add_argument('--tolerance', help='Deviation from expected length allowed', type = int, default = 10)
    parser.add_argument('--data_dir', help='Location of the data directory of package', type = str, default="./data")
    parser.add_argument('--output', required=True, help='Prefix of output files  (required)', type = str)

    args = parser.parse_args()
    
    return args

def search_for_gene(ref_in,name,gene_length,tol,correct_length,gene_rower):

    if not correct_length:
        gene_row = ref_in['attributes'].str.contains(name)
        gene_row_indy = gene_row.where(gene_row == True)
        gene_row_indy = gene_row_indy.index[gene_row_indy == True].tolist()
        gene_rower = ref_gff_tsv.iloc[gene_row_indy]
        if gene_rower.empty == False:
            gene_len = [abs(int(gene_rower.iloc[0,4]) - int(gene_rower.iloc[0,3]))]

            overhang = [gene_length - tol, gene_length + tol]

            if overhang[0] <= gene_len[0] <= overhang[1]:
                correct_length = True
            else:
                sys.stderr.write('Found gene' + name + ' but wrong length: ' + str(gene_len[0]) + ', expected: ' + str(gene_length) + '\n')

    return correct_length,gene_rower

def get_aln_pos_from_ref(hmm_aln,pos,offset):
    
    pos = pos - offset
    ref_pos = pos + 1
    upstream_length = 0
    while upstream_length < pos:
        upstream_frag = hmm_aln[0,0:ref_pos].seq
        upstream_length = len(upstream_frag) - upstream_frag.count('.')
        ref_pos = ref_pos + 1
    return(ref_pos - 1)

def hmm_search_for_gene(fasta,gene, aa_dir_name, data_dir):
    
    # data structure
    readingframes = []
    
    # from https://www.biostars.org/p/183616/
    # load sequence

    tic_aa_creator = time.perf_counter()

    for record in SeqIO.parse(fasta, "fasta"):


        # process file name expecting no extra .
        gff_base = os.path.basename(fasta)
        gff_base = re.sub("\..*[a-z,A-Z,0-9].*$","",gff_base)
        aa_base_name = aa_dir_name + "/" + gff_base

        # Create three reading frames in forward direction, offset 0, 1, 2
        record_rc = record.reverse_complement()
        readingframes.extend([Seq.translate(record.seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in range(3)])
        readingframes.extend([Seq.translate(record_rc.seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in range(3)])

        results = []
        for frame in readingframes:
            for peptide in frame.split('*'): #Split translation over stopcodons
                if len(peptide) > 30:
                    results.append(peptide)

        #Write length and translation to file
        #Use PotentialORFs.txt as output, can be changed
        with open(aa_base_name + '.aa', 'w') as output:
            for n,peptide in enumerate(results):
                output.write(">{}\n{}\n".format('peptide_' + str(n), peptide))

    toc_aa_creator = time.perf_counter()

    tic_hmm_run = time.perf_counter()
    # run HMM
    hmm_output_fn = aa_base_name + '.' + gene + '.hsp.out'
    hmm_proc_out = 1
    while hmm_proc_out != 0:
        hmm_proc_out = subprocess.check_call('hmmsearch ' + data_dir + '/' + gene + '.hmm ' + aa_base_name + '.aa > ' + hmm_output_fn, shell = True)

    # parse HMM
    hmm_output = list(SearchIO.parse(hmm_output_fn, 'hmmer3-text'))
    hmm_best_hsp_index = None
    hmm_best_hsp_bitscore = 0
    for i,hsp in enumerate(hmm_output[0]):
        if hsp.bitscore > hmm_best_hsp_bitscore:
            hmm_best_hsp_bitscore = hsp.bitscore
            hmm_best_hsp_index = i

    if hmm_best_hsp_index is None:
        return "Missing_HMM", gff_base


    hmm_best_hsp = hmm_output[0][hmm_best_hsp_index][0]
    hmm_offset = hmm_best_hsp.query_start
    hmm_aln = hmm_best_hsp.aln
    subprocess.call('rm ' + hmm_output_fn, shell = True)

    toc_hmm_run = time.perf_counter()

    if gene == 'folP':
        aln_start = get_aln_pos_from_ref(hmm_aln,57,hmm_offset)
        aln_end = get_aln_pos_from_ref(hmm_aln,67,hmm_offset)
        hmm_match = hmm_aln[0,aln_start:aln_end].seq
        query_match = hmm_aln[1,aln_start:aln_end].seq
        gap_count = hmm_match.count('.') - query_match.count('.')
        
        if gap_count > 0:
            print(gff_base + '\tSulphamethoxazole resistant')
            status = "R"
        else:
            print(gff_base + '\tSulphamethoxazole sensitive')
            status = "S"
    
    elif gene == 'dhfR':
        aln_start = get_aln_pos_from_ref(hmm_aln,98,hmm_offset)
        hmm_match = hmm_aln[0,aln_start:(aln_start+1)].seq
        query_match = hmm_aln[1,aln_start:(aln_start+1)].seq
        
        if query_match.upper() == "L":
            print(gff_base + '\tTrimethoprim resistant (' + query_match.upper() + ')')
            status = "R"
        elif query_match.upper() == "I":
            print(gff_base + '\tTrimethoprim sensitive (' + query_match.upper() + ')')
            status = "S"
        else:
            print(gff_base + '\tTrimethoprim unknown: ' + str(query_match).upper())
            status = "NA"

    print("Took this long for ORF finder: %s" % (toc_aa_creator - tic_aa_creator))
    print("Took this long for HMM run: %s" % (toc_hmm_run - tic_hmm_run))

    return status, gff_base

def extract_fasta_from_gff(gff_fn):

    # check file name
    fa_fn = gff_fn
    if gff_fn.endswith('.gff'):
        fa_fn = fa_fn.replace('.gff','.fa')
    else:
        sys.stderr.write('Does not appear to be a GFF file: ' + gff_fn + '\n')
        exit(1)

    # write output
    in_fa_region = False
    with open(gff_fn,'r') as gff_file, open(fa_fn,'w') as fa_file:
        for line in gff_file.readlines():
            if line.startswith('>'):
                in_fa_region = True
            if in_fa_region:
                fa_file.write(line)
    
    # return name
    return fa_fn

def process_pbp(isolate,protein_fasta,tpd_start,tpd_end,tpd_lab,results_csv,k):
    protein_rec = list(SeqIO.parse(protein_fasta, format="fasta"))
    prot = protein_rec[0].seq

    tpd_string = str(prot[tpd_start: tpd_end])

    tpd_string = Seq(tpd_string)#, generic_protein)

    prot_tpd_id = isolate + "_" + tpd_lab + "_TPD"
    tpd_protein_string = SeqIO.SeqRecord(tpd_string, id=prot_tpd_id)


    prot_file = isolate + "_" + tpd_lab  + ".prot"
    with open(prot_file, "w+") as output_handle:
        SeqIO.write(tpd_protein_string, output_handle, "fasta")

    top_id = str(results_csv.iloc[0,1])
    top_id = top_id[2:]

    

    if bassio_nameo != "22841_3#15.contigs_velvet.fa.gff":

        rm_command = "rm " + protein_fasta + " " + protein_csv
        os.system(rm_command)
    else:
        print(sstart, send)
        print(tpd_start, tpd_end)

    print("Generating CSV: %s%%" % round((k / len(gff_lines) * 100)))

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

    gff_files = open(files_for_input.gff, "r")
    gff_lines = gff_files.read().splitlines()
    
    fasta_lines = []
    if files_for_input.fasta is None:
        for gff in gff_lines:
            fasta_lines.append(extract_fasta_from_gff(gff))
    else:
        fasta_lines = open(files_for_input.fasta,"r")
        fasta_lines = fasta_lines.read().splitlines()

    gene_name = files_for_input.pbp
    input_gene_names = files_for_input.gene.split(',') # allow a list of alternative names
    input_gene_names.append(gene_name)
    input_gene_names.append(gene_name.lower())
    all_gene_names = set(input_gene_names)
    gene_length = files_for_input.tlength

    if files_for_input.pbp not in ['pbp1a','pbp2b','pbp2x']:
        aa_dir_name = "./" + re.sub(".csv","",files_for_input.output) + "_aa_dir"
        if not os.path.exists(aa_dir_name):
            os.mkdir(aa_dir_name)
        all_gene_names = [files_for_input.pbp]

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
    skip = False
    toc_setup = time.perf_counter()
    print("Took this long for initial set up: %s" % (toc_setup - tic_setup))
    print("Beginning iso run")
    for k,(gff_file,fasta_file) in enumerate(zip(gff_lines,fasta_lines)):
        print("On isolate %s of %s" % (k, len(gff_lines)))
        tic_iso_run = time.perf_counter()
        bassio_nameo = os.path.basename(gff_file)
        
        ref_gff_tsv = pandas.read_csv(gff_file, sep='\t',
                                      names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                             'attributes'],
                                      header=None)
        
        if ref_gff_tsv['attributes'].isnull().all():
            missing_isolates.append(bassio_nameo)
            print(("This isn't a recognised GFF: " + bassio_nameo))
            continue
        
        gene_rower = None
        correct_length = False
        for gene in all_gene_names:
            if not correct_length:
                #print('Searching for name ' + gene)
                if files_for_input.pbp in ['pbp1a','pbp2b','pbp2x']:
                    correct_length,gene_rower = search_for_gene(ref_gff_tsv,
                                                                gene,
                                                                gene_length,
                                                                files_for_input.tolerance,
                                                                correct_length,
                                                                gene_rower)
                else:

                    current_res, current_isolate = hmm_search_for_gene(fasta_file, files_for_input.pbp, aa_dir_name, files_for_input.data_dir)
                    gene_id.append(current_res)
                    isolate_name.append(current_isolate)
                    skip = True
                    toc_iso_run = time.perf_counter()
                    print("Took this long for isolate %s, overall: %s" %(current_isolate, (toc_iso_run - tic_iso_run)))
                    continue
                if correct_length:
                    print('Found gene ' + gene + ' in ' + bassio_nameo + '\n')

        if skip:
            continue

        if gene_rower.empty:
            
            print("No gene hit in this file:", bassio_nameo)
            missing_isolates.append(bassio_nameo)
            
        else:

            # extract hit information
            contig_id = gene_rower.iloc[0,0]
            gene_start = int(gene_rower.iloc[0,3])
            gene_end = int(gene_rower.iloc[0,4])
            strand = str(gene_rower.iloc[0, 6])

            # identify contigs within assembly
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id == contig_id:

                    correct_contig = record.seq

            
            # identify gene sequence within contig
            gene_string = str(correct_contig[(gene_start - 1):gene_end])
            if strand == "-":

                reverso = Seq(gene_string)#, generic_dna)

                gene_string = str(reverso.reverse_complement())

            isolate = re.split("\.", bassio_nameo)[0]

            # translate gene to protein
            

            protein_string = Seq(gene_string).translate()#, generic_dna).translate()
            protein_string = SeqIO.SeqRecord(protein_string, id = bassio_nameo)
            
            
            print(bassio_nameo + ".fasta")
            # write protein sequence out for BLAST
            with open((bassio_nameo + files_for_input.pbp + ".fasta"), "w+") as output_handle:
                SeqIO.write(protein_string, output_handle, "fasta")

            #######################################################################
            ## Now we've got our protein string we'll use blast to search for the #
            ## closest hit in the cdc alignment file ##############################
            #######################################################################

            protein_fasta = bassio_nameo + files_for_input.pbp + ".fasta"
            protein_csv = bassio_nameo + files_for_input.pbp  + ".csv"

            align_path = files_for_input.aln
            basename_aligno = os.path.basename(files_for_input.aln)
            basename_blastdb = basename_aligno + ".phr"

            if not os.path.isfile(basename_blastdb):
                db_command = "makeblastdb -in " + align_path + " -dbtype prot -out " + basename_aligno
                subprocess.call(db_command, shell = True)
            
            blasty_cmd = "blastp -query " + protein_fasta + " -db " + basename_aligno + " -out " + protein_csv + " -outfmt 10"
            subprocess.call(blasty_cmd, shell = True)

            #######################################################################
            ## Now we've got the blast results, lets extract the top hit and then #
            ## assign this isolates name as that particular allele type ###########
            #######################################################################

            results_csv = pandas.read_csv(protein_csv,
                                          names=['qid', 'sid', 'pid', 'align',
                                                 'gap', 'mismatch', 'qstart', 'qend',
                                                 'sstart','send','eval','bitscore'])


            results_csv = results_csv.sort_values(['bitscore'], ascending = False)
            if results_csv.empty:
                if bassio_nameo == "14673_2#1.contigs_velvet.fa.gff":
                    results_csv.to_csv("14673_2#1.contigs_velvet.fa.gff.csv",
                                       index=False)
                print("No Blast results for this isolate:" , bassio_nameo)
                missing_isolates.append(bassio_nameo)
            else:
                top_res = results_csv.iloc[0]
                sstart = top_res.iloc[6] - 1
                send = top_res.iloc[7]

                if files_for_input.pbp == "pbp1a":
                    tpd_start = sstart
                    tpd_end = sstart + 252
                    tpd_lab = "pbp1a"
                elif files_for_input.pbp == "pbp2b":
                    tpd_start = sstart
                    tpd_end = sstart + 277
                    tpd_lab = "pbp2b"
                elif files_for_input.pbp == "pbp2x":
                    tpd_start = sstart + 60
                    tpd_end = tpd_start + 299
                    tpd_lab = "pbp2x"
                    if isolate == "11511_7#11":
                        print(sstart, tpd_start,"Now for the ends \n", send, tpd_end)
                        print(top_res)
                else:
                    sys.exit("For TPD extraction please use one of the three following gene names: pbp1a, pbp2b or pbp2x")

                isolate_name.append(isolate)
                top_pbp_hit = process_pbp(isolate,protein_fasta,tpd_start,tpd_end,tpd_lab,results_csv,k)
                gene_id.append(top_pbp_hit)

    out_df = pandas.DataFrame()
    out_df['isolate_id'] = pandas.Series(isolate_name)
    if files_for_input.pbp in ['pbp1a','pbp2b','pbp2x']:
        out_df['top_id'] = pandas.Series(gene_id, index=out_df.index)
    else:
        out_df['Resistance'] = pandas.Series(gene_id, index=out_df.index)
        subprocess.call("rm -r *_aa_dir", shell=True)

    # rm_command = "rm " + protein_fasta + " " + protein_csv
    # os.system(rm_command)
    
    if len(missing_isolates) > 0:
        with open(("./missing_" + files_for_input.pbp + "_ids.txt"), mode='wt', encoding='utf-8') as myfile:
            myfile.write('\n'.join(missing_isolates) + '\n')
    out_df.to_csv(files_for_input.output,
                  index=False)













