#!/bin/bash

## Shell script to extract pbp genes from a list of gffs and fastas

set -e
echo "This is the gff list $1"
echo "This is the fasta list $2"
echo "This is the out csv name $3"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
parentdir="$(dirname "$DIR")"
pythondir="${parentdir}/python/"
rdir="${parentdir}/R/"
datadir="${parentdir}/data/"
echo $-
if [[ $- == *i* ]]
then
  echo "Interactive"
  read -p "Are you sure (Y/N)? " -n 1 -r
  echo    # (optional) move to a new line
  if [[ $REPLY =~ ^[Nn]$ ]]
  then
      echo "exiting script now"
    else
    ## Lets run through the pbp gene getting pbp1a, 2b 2x
    echo "Extracting pbp1a/ponA now"
    python "${pythondir}pen_checker_cdc.py" \
    --gff $1 --fasta $2 --pbp pbp1a --gene ponA,pbp1A --aln "${datadir}pbp1a.faa" --tlength 2159 --output pbp1a_cdc_mic_type.csv --tolerance 100

    echo "Extracting pbp2b/penA now"
    python "${pythondir}pen_checker_cdc.py" \
    --gff $1 --fasta $2 --pbp pbp2b --gene penA,pbp2B --aln "${datadir}pbp2b.faa" --tlength 2042 --output pbp2b_cdc_mic_type.csv --tolerance 100

    echo "Extracting pbp2x/pbpX now"
    python "${pythondir}pen_checker_cdc.py" \
    --gff $1 --fasta $2 --pbp pbp2x --gene pbpX,pbp2X --aln "${datadir}pbp2x.faa" --tlength 2252 --output pbp2x_cdc_mic_type.csv --tolerance 100

    ## Create the lists of prots for the Rscript to work on

    ls -d $PWD/*pbp1a.prot > pbp1a_prot_list
    ls -d $PWD/*pbp2b.prot > pbp2b_prot_list
    ls -d $PWD/*pbp2x.prot > pbp2x_prot_list

    ## Lets run the Rscript to create the AA df

    Rscript --vanilla "${rdir}AA_df_creator.R" \
    pbp1a_prot_list pbp2b_prot_list pbp2x_prot_list aa_df.csv

    ## Lets run the Rscript to predict the cats for the model
    Rscript --vanilla "${rdir}RF_run.R" \
    aa_df.csv "${datadir}cdc_seqs_df.csv" 3 $3

  fi
else
  echo "non interactive"
  ## Lets run through the pbp gene getting pbp1a, 2b 2x
    echo "Extracting pbp1a/ponA now"
    python "${pythondir}pen_checker_cdc.py" \
    --gff $1 --pbp pbp1a --gene ponA,pbp1A,"pbp 1 A","pbp1 A","pbp 1 a","pbp1 a" --aln "${datadir}pbp1a.faa" --tlength 2159 --output pbp1a_cdc_mic_type.csv --tolerance 500

    echo "Extracting pbp2b/penA now"
    python "${pythondir}pen_checker_cdc.py" \
    --gff $1 --pbp pbp2b --gene penA,pbp2B,"pbp 2 B","pbp2 B","pbp 2 b","pbp2 b" --aln "${datadir}pbp2b.faa" --tlength 2042 --output pbp2b_cdc_mic_type.csv --tolerance 650

    echo "Extracting pbp2x/pbpX now"
    python "${pythondir}pen_checker_cdc.py" \
    --gff $1 --pbp pbp2x --gene pbpX,pbp2X,"pbp 2 X","pbp2 X","pbp 2 x","pbp2 x","pbp2x" --aln "${datadir}pbp2x.faa" --tlength 2252 --output pbp2x_cdc_mic_type.csv --tolerance 500

    ## Create the lists of prots for the Rscript to work on

    ls -d $PWD/*pbp1a.prot > pbp1a_prot_list
    ls -d $PWD/*pbp2b.prot > pbp2b_prot_list
    ls -d $PWD/*pbp2x.prot > pbp2x_prot_list

    ## Lets run the Rscript to create the AA df

    Rscript --vanilla "${rdir}AA_df_creator.R" \
    pbp1a_prot_list pbp2b_prot_list pbp2x_prot_list aa_df.csv

    ## Lets run the Rscript to predict the cats for the model
    Rscript --vanilla "${rdir}RF_run.R" \
    aa_df.csv "${datadir}cdc_seqs_df.csv" 3 $3

    rm *.prot
    rm pbp1a.faa.p* pbp2b.faa.p* pbp2x.faa.p*

fi


