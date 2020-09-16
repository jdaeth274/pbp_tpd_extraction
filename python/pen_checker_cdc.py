from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
import pandas
import numpy
import re
import os
import sys
import argparse
import subprocess

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
    parser.add_argument('--fasta', required=True, help='List of FASTA files (required)', type = str)
    parser.add_argument('--pbp', required=True, help='Specify PBP being analysed  (required)', type = str, choices=['pbp1a', 'pbp2b', 'pbp2x'])
    parser.add_argument('--gene', required=True, help='Comma-separated list of gene synonyms (required)', type = str)
    parser.add_argument('--aln', required=True, help='Gene alignment  (required)', type = str)
    parser.add_argument('--tlength', required=True, help='Target length of TPD domain  (required)', type = int)
    parser.add_argument('--tolerance', help='Deviation from expected length allowed', type = int, default = 10)
    parser.add_argument('--output', required=True, help='Prefix of output files  (required)', type = str)

    args = parser.parse_args()
    
    return args

def search_for_gene(ref_in,name,gene_length,tol,correct_length):

    if not correct_length:
        gene_row = ref_in['attributes'].str.contains(name)
        gene_row_indy = gene_row.where(gene_row == True)
        gene_row_indy = gene_row_indy.index[gene_row_indy == True].tolist()
        ## If multiple paralogs with _ this just assumes the first found is the correct gene.
        if len(gene_row_indy) > 1:
            gene_row_indy = [gene_row_indy[0]]
        if not gene_row_indy == False:
            gene_rower = ref_in.iloc[gene_row_indy]
        if gene_rower.empty == False:
            if isinstance(gene_rower, pandas.Series):
                gene_len = [abs(int(gene_rower.iloc[4]) - int(gene_rower.iloc[3]))]
                gene_rower = gene_rower.to_frame().reset_index()
            else:
                gene_len = [abs(int(gene_rower.iloc[0,4]) - int(gene_rower.iloc[0,3]))]

            overhang = [gene_length - tol, gene_length + tol]

            if overhang[0] <= gene_len[0] <= overhang[1]:
                correct_length = True
            else:
                sys.stderr.write('Found gene ' + name + ' but wrong length: ' + str(gene_len[0]) + ', expected: ' + str(gene_length) + '\n')

    return correct_length,gene_rower

#########
# BEGIN #
#########

if __name__ == '__main__':

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

    fasta_lines = open(files_for_input.fasta,"r")
    fasta_lines = fasta_lines.read().splitlines()



    gene_name = files_for_input.pbp
    input_gene_names = files_for_input.gene.split(',') # allow a list of alternative names
    input_gene_names.append(gene_name)
    input_gene_names.append(gene_name.lower())
    all_gene_names = set(input_gene_names)
    gene_length = files_for_input.tlength
    
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

    for k in range(len(gff_lines)):

        bassio_nameo = os.path.basename(gff_lines[k])
        
        ref_gff_tsv = pandas.read_csv(gff_lines[k], sep='\t',
                                      names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                             'attributes'],
                                      header=None)
        gene_rower = None
        correct_length = False
        for gene in all_gene_names:
            if not correct_length:
                #print('Searching for name ' + gene)
                correct_length,gene_rower = search_for_gene(ref_gff_tsv,gene,gene_length,files_for_input.tolerance,correct_length)
                #if correct_length:
                #    sys.stderr.write('Found gene ' + gene + ' in ' + bassio_nameo + '\n')

        if gene_rower.empty:
            
            print("No gene hit in this file:", bassio_nameo)
        else:
            bassio_nameo = os.path.basename(gff_lines[k])
            contig_num = gene_rower.iloc[0,0]

            gene_start = int(gene_rower.iloc[0,3])
            gene_end = int(gene_rower.iloc[0,4])
            strand = str(gene_rower.iloc[0, 6])

            #print(contig_num)

            contig_num = re.findall(r'\d+', contig_num)[-1]
            #print(contig_num)
            #contig_num = re.split("tig", contig_num)[-1]
            contig_num = int(contig_num)



            fasta_file = fasta_file_searcher(gff_lines[k], fasta_lines)



            records = list(SeqIO.parse(fasta_file, "fasta", alphabet=generic_dna))
            if bassio_nameo == "INV200.contigs_velvet.fa.gff":
                correct_contig = records[0].seq
            else:
                correct_contig = records[contig_num -1].seq

            gene_string = str(correct_contig[(gene_start - 1):gene_end])
            if strand == "-":
                reverso = Seq(gene_string, generic_dna)
                gene_string = str(reverso.reverse_complement())

            isolate = re.split("\.", bassio_nameo)[0]

            protein_string = Seq(gene_string, generic_dna).translate()
            protein_string = SeqIO.SeqRecord(protein_string, id = bassio_nameo)



            with open((bassio_nameo + ".fasta"), "w+") as output_handle:
                SeqIO.write(protein_string, output_handle, "fasta")

            if bassio_nameo == "14673_2#1.contigs_velvet.fa.gff":
                with open("14673_2#1.contigs_velvet.fa.gff,fasta", "w+") as output_handle:
                    SeqIO.write(protein_string, output_handle, "fasta")

            #######################################################################
            ## Now we've got our protein string we'll use blast to search for the #
            ## closest hit in the cdc alignment file ##############################
            #######################################################################

            protein_fasta = bassio_nameo + ".fasta"
            protein_csv = bassio_nameo + ".csv"

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

                protein_rec = list(SeqIO.parse(protein_fasta, format="fasta",
                                               alphabet=generic_protein))
                prot = protein_rec[0].seq

                tpd_string = str(prot[tpd_start: tpd_end])
                tpd_string = Seq(tpd_string, generic_protein)
                prot_tpd_id = isolate + "_" + tpd_lab + "_TPD"
                tpd_protein_string = SeqIO.SeqRecord(tpd_string, id=prot_tpd_id)


                prot_file = isolate + "_" + tpd_lab  + ".prot"
                with open(prot_file, "w+") as output_handle:
                    SeqIO.write(tpd_protein_string, output_handle, "fasta")

                top_id = str(results_csv.iloc[0,1])
                top_id = top_id[2:]
                isolate_name.append(isolate)
                gene_id.append(top_id)

                if bassio_nameo != "22841_3#15.contigs_velvet.fa.gff":

                    rm_command = "rm " + protein_fasta + " " + protein_csv
                    os.system(rm_command)
                else:
                    print(sstart, send)
                    print(tpd_start, tpd_end)

                print("Generating CSV: %s%%" % round((k / len(gff_lines) * 100)))


    out_df = pandas.DataFrame()
    out_df['isolate_id'] = pandas.Series(isolate_name)
    out_df['top_id'] = pandas.Series(gene_id, index=out_df.index)

    # rm_command = "rm " + protein_fasta + " " + protein_csv
    # os.system(rm_command)


    out_df.to_csv(files_for_input.output,
                  index=False)













