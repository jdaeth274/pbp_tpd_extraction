#!/bin/bash

## Shell script to extract pbp genes from a list of gffs and fastas

set -e

if [[ $# == 0 ]]
then
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  echo "Missing files check the usage below"
  echo "Requires: [-f|--fasta-list] [-g|--gene (pbp/dhfr/folp)] [-w|--whole-gene (Y/N)] [-o|--out-prefix] [-t|--threads]"
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  exit
fi

POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
  -f|--fasta-list)
  FASTA="$2"
  shift
  shift
  ;;
  -g|--gene)
  GENE="$2"
  shift
  shift
  ;;
  -w|--whole-gene)
  WHOLE="$2"
  shift
  shift
  ;;
  -o|--out-prefix)
  OUT="$2"
  shift
  shift
  ;;
  -t|--threads)
  THREADS="$2"
  shift
  shift
  ;;
esac
done

set -- "${POSITIONAL[@]}"

if [ $GENE != "pbp" ] && [ $GENE != "dhfr" ] && [ $GENE != "folp" ]
then
  echo "Gene needs to be one of the following three: "
  echo "pbp, dhfr or folp"
  echo "exiting script"
  exit
fi

echo "This is the fasta list: ${FASTA}"
echo "This is the gene(s) to search: ${GENE}"
echo "Perform whole gene extraction? ${WHOLE}"
echo "This is the out prefix: ${OUT}"
echo "This is the number of threads to use: ${THREADS}"
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

    if [ $WHOLE == "Y" ] || [ $WHOLE == "y" ] || [ $WHOLE == "yes" ] || [ $WHOLE == "YES" ] || [ $WHOLE == "Yes" ]
    then
      echo "Extracting whole sequences of pbp genes now"
      python "${pythondir}gene_extraction.py" \
      --seqs $FASTA --data_dir $datadir --output $OUT --cores $THREADS

    elif [ $GENE == "pbp" ]
    then
      echo "Extracting pbp1a now"
      python "${pythondir}orfipy_test.py" \
      --seqs $FASTA --gene pbp1a --tlength 720 --aln "${datadir}pbp1a.faa" --tolerance 100 \
      --data_dir $datadir --output "${OUT}_1a_isos.txt" --cores $THREADS

      echo "Extracting pbp2b now"
      python "${pythondir}orfipy_test.py" \
      --seqs $FASTA --gene pbp2b --tlength 681 --aln "${datadir}pbp2b.faa" --tolerance 100 \
      --data_dir $datadir --output "${OUT}_2b_isos.txt" --cores $THREADS

      echo "Extracting pbp2x now"
      python "${pythondir}orfipy_test.py" \
      --seqs $FASTA --gene pbp2x --tlength 681 --aln "${datadir}pbp2x.faa" --tolerance 100 \
      --data_dir $datadir --output "${OUT}_2x_isos.txt" --cores $THREADS

      ls -d $PWD/*pbp1a.prot > pbp1a_prot_list
      ls -d $PWD/*pbp2b.prot > pbp2b_prot_list
      ls -d $PWD/*pbp2x.prot > pbp2x_prot_list

      ## Lets run the Rscript to create the AA df

      Rscript --vanilla "${rdir}AA_df_creator.R" \
      pbp1a_prot_list pbp2b_prot_list pbp2x_prot_list aa_df.csv

      ## Lets run the Rscript to predict the cats for the model
      Rscript --vanilla "${rdir}RF_run.R" \
      aa_df.csv "${datadir}cdc_seqs_df.csv" 3 "${OUT}_res_preds.csv"

    elif [ $GENE == "dhfr" ]
    then
      echo "Checking Trimethoprim (dhfR) resistance"
      python "${pythondir}orfipy_test.py" \
      --seqs $FASTA --gene dhfR --tlength 2159 --output "${OUT}.csv" --tolerance 100 \
      --data_dir $datadir --cores $THREADS

    elif [ $GENE == "folp" ]
    then
      echo "Checking sulfamethoxazole (folP) resistance"
      python "${pythondir}orfipy_test.py" \
      --seqs $FASTA --gene folP --tlength 2159 --output "${OUT}.csv" --tolerance 100 \
      --data_dir $datadir --cores $THREADS
    else
      echo "Unrecognised command please check usage:"
      echo ""
      echo "bash pbp_gene_extraction.sh [-f|--fasta-list] [-g|--gene (pbp/dhfr/folp)] [-w|--whole-gene (Y/N)] [-o|--out-prefix] [-t|--threads]"
      echo ""
      exit

    fi

    if [ -d "${OUT}_prots" ]
    then
      rm -r "${OUT}_prots"
      mkdir "${OUT}_prots"
    else
      mkdir "${OUT}_prots"
    fi

    mv *.prot "${OUT}_prots"



  fi
else
  echo "non interactive"
  ## Lets run through the pbp gene getting pbp1a, 2b 2x
    if [ $WHOLE == "Y" ] || [ $WHOLE == "y" ] || [ $WHOLE == "yes" ] || [ $WHOLE == "YES" ] || [ $WHOLE == "Yes" ]
    then
      echo "Extracting whole sequences of pbp genes now"
      python "${pythondir}gene_extraction.py" \
      --seqs $FASTA --data_dir $datadir --output $OUT --cores $THREADS

    elif [ $GENE == "pbp" ]
    then
      echo "Extracting pbp1a now"
      python "${pythondir}orfipy_test.py" \
      --seqs $FASTA --gene pbp1a --tlength 720 --aln "${datadir}pbp1a.faa" --tolerance 100 \
      --data_dir $datadir --output "${OUT}_1a_isos.txt" --cores $THREADS

      echo "Extracting pbp2b now"
      python "${pythondir}orfipy_test.py" \
      --seqs $FASTA --gene pbp2b --tlength 681 --aln "${datadir}pbp2b.faa" --tolerance 100 \
      --data_dir $datadir --output "${OUT}_2b_isos.txt" --cores $THREADS

      echo "Extracting pbp2x now"
      python "${pythondir}orfipy_test.py" \
      --seqs $FASTA --gene pbp2x --tlength 681 --aln "${datadir}pbp2x.faa" --tolerance 100 \
      --data_dir $datadir --output "${OUT}_2x_isos.txt" --cores $THREADS

      ls -d $PWD/*pbp1a.prot > pbp1a_prot_list
      ls -d $PWD/*pbp2b.prot > pbp2b_prot_list
      ls -d $PWD/*pbp2x.prot > pbp2x_prot_list

      ## Lets run the Rscript to create the AA df

      Rscript --vanilla "${rdir}AA_df_creator.R" \
      pbp1a_prot_list pbp2b_prot_list pbp2x_prot_list aa_df.csv

      ## Lets run the Rscript to predict the cats for the model
      Rscript --vanilla "${rdir}RF_run.R" \
      aa_df.csv "${datadir}cdc_seqs_df.csv" 3 "${OUT}_res_preds.csv"

    elif [ $GENE == "dhfr" ]
    then
      echo "Checking Trimethoprim (dhfR) resistance"
      python "${pythondir}orfipy_test.py" \
      --seqs $FASTA --gene dhfR --tlength 2159 --output "${OUT}.csv" --tolerance 100 \
      --data_dir $datadir --cores $THREADS

    elif [ $GENE == "folp" ]
    then
      echo "Checking sulfamethoxazole (folP) resistance"
      python "${pythondir}orfipy_test.py" \
      --seqs $FASTA --gene folP --tlength 2159 --output "${OUT}.csv" --tolerance 100 \
      --data_dir $datadir --cores $THREADS
    else
      echo "Unrecognised command please check usage:"
      echo ""
      echo "bash pbp_gene_extraction.sh [-f|--fasta-list] [-g|--gene (pbp/dhfr/folp)] [-w|--whole-gene (Y/N)] [-o|--out-prefix] [-t|--threads]"
      echo ""
      exit

    fi

    if [ -d "${OUT}_prots" ]
    then
      rm -r "${OUT}_prots"
      mkdir "${OUT}_prots"
    else
      mkdir "${OUT}_prots"
    fi

    mv *.prot "${OUT}_prots"
    #rm *.prot
    #rm pbp1a.faa.p* pbp2b.faa.p* pbp2x.faa.p*

fi


