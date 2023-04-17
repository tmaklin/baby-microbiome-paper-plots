#!/bin/bash

## Extract the lineages identified at each time point and by individuali

# NOTE:
## Should be run from the source directory with `Scripts/extract_time_points.sh`

set -euo pipefail

## Remove the header lines
sed '1d' wgs_demix_check_high_confidence.tsv | sort -k1 > tmp1.txt
sed '1d' wgs_meta_delivery.tsv | sort -k1 > tmp2.txt

## Add header for the extracted columns
touch tmp.tsv
echo -e "accession\tdelivery_mode\tindividual\ttime_point\tlineage\trelative_abundance\tdemix_check_score" > tmp.tsv

## Extract values
join -1 1 -2 1 tmp2.txt tmp1.txt | tr ' ' '\t' | cut -f1,2,3,4,5,6,8 >> tmp.tsv

## Rename the PopX clusters to (Species)_SC_ST format
for f in poppunk_cluster_info/*.tsv; do while read line; do cluster=$(echo $line | cut -f1 -d' '); newcluster=$(echo $line | cut -f2 -d' '); sed -i "s/$cluster[[:space:]]/$newcluster\t/g" tmp.tsv; done < <(cat $f | sed '1d'); done | head

mv tmp.tsv babybiome_lineages_by_time_point.tsv

