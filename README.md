## Extraction the pbp TPD domains from Strep pneumo gffs ##

# Installation #

This requires conda, please install conda first [here](docs.conda.io/projects/conda/en/latest/user-guide/install)
Once installed clone the repo
`git clone https://github.com/jdaeth274/pbp_tpd_extraction`

Then use the environment.yml file to install the dependencies with conda

`conda env create --file=environment.yml`

Activate this environment using 

`conda activate pbp_tpd_env`

# Usage # 

For pbp extraction, once in the pbp_tpd_env created from above use the following script command
to run and extract a csv of the isolate name and the pbp resistance category. 

`bash ./bash/pbp_gene_extraction.sh gff_list.txt fasta_list.txt out_csv_name.csv`

Here the gff_list is a txt file with the location of a gff file on each line.
The fasta list is complementary to this, with the same fasta file on each line, 
as corresponds to the gff list.


