#!/bin/bash
## Script for creating a multifasta file from assemblies to use with
## Themisto/mSWEEP/mGEMS.
##
## Usage:
##   create_reference.sh assembly_info.tsv gap_length > output.fasta
## Input:
##   assembly_info.tsv: tab-separated file containing the names, cluster assignments,
##                      and paths of each assembly to include in the multifasta. The
##                      file should be in the following format (without the leading #'s):
##                      #id        cluster    assembly
##                      #seq_1     clu_1      /path/to/assembly/seq_1.fasta.gz
##                      #seq_2     clu_1      /path/to/assembly/seq_2.fasta
##   gap_length: if the assemblies are fragmented, this value tells how long gap will be
##               placed between the fragments. Recommended: 2x read length.
## Output:
##   (plaintext) multifasta file will be printed by the program;
##   usage example saves the output to the file `output.fasta`.
##
set -euxo pipefail

gap=''
for (( i = 0; i < $2; i++ )); do
    gap=$gap'-'
done

datafile=tmp_$RANDOM""_information.tsv
sed '1d' $1 > $datafile

while read line; do
    id=$(echo $line | cut -f1 -d' ')
    assembly=$(echo $line | cut -f3 -d' ')
    echo ">"$id ## new sequences in multifasta start with `>`
    ## Print the assembly if its gzipped (--force to allow printing plaintext files)
    ## subcommands: sed "s/^>.*$/$gap/g" - replace new contigs starting with `>` with a gap denoted by -'s.
    ##              sed '1d' - removes the leading gap.
    ##              tr -d '\n' - remove line breaks.
    ##              fold -c79 - break the sequence in lines that contain 80 characters each
    ##              sed '$s/$/\n/' - add a newline at the end of the sequence.
    zcat --force $assembly | sed "s/^>.*$/$gap/g" | sed '1d' | tr -d '\n' | fold -c79 | sed '$s/$/\n/'
done < $datafile

rm $datafile
