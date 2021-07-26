# Extraction of the pbp TPD domains from Strep pneumo gffs #

## Installation ##

This requires conda, please install conda first [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
Once installed clone the repo:

    git clone https://github.com/jdaeth274/pbp_tpd_extraction

Then use the environment.yml file to install the dependencies with conda

    conda env create --file=environment.yml

Activate this environment using 

    conda activate pbp_tpd_env

## Usage ##

### Penicillin resistance ###
To extract the pbp genes (_pbp1a, pbp2b, pbp2x_) Transpeptidase domains (TPDs) and run the RF prediction model use the 
following command: 

    bash pbp_gene_extraction.sh --fasta-list seqeunce_list.txt --gene pbp --whole-gene N --out-prefix test_out --threads 1

This will produce a directory containing the extracted proteins, a csv file `aa_df.csv` of the extracted TPD domains 
per isolate for the three genes, a prediction csv and separate files containing the isolates with missing domains for 
each gene. 

To extract the whole gene seqeunces for the pbp genes only use:

    bash pbp_gene_extraction.sh --fasta-list seqeunce_list.txt --gene pbp --whole-gene Y --out-prefix test_out --threads 1

### Co-trimoxazole resistance ###

To extract the resistance profiles for trimethoprim use the following command: 

    bash pbp_gene_extraction.sh --fasta-list seqeunce_list.txt --gene dhfr --whole-gene N --out-prefix test_out --threads 1

To extract the resistance profiles for sulfamethoxazole use the following command:

    bash pbp_gene_extraction.sh --fasta-list seqeunce_list.txt --gene folp --whole-gene N --out-prefix test_out --threads 1







