#!/bin/bash

inFile=/Users/shraddhapai/Documents/Research/Epigenetics/NARSAD2014/anno/gencode.v27lift37.annotation.gtf

echo "getting genes"
cat $inFile | awk '{if ($3=="gene"){print}}' > ${inFile}.genes
echo "getting coords"
cat ${inFile}.genes | cut -f 1,4,5,7 > ${inFile}.chroms
echo "getting names"
cat ${inFile}.genes | cut -f 9 |  awk '{gsub(/["]/,"");}1' > ${inFile}.geneids #head | perl -pe '/gene_id ([^"]+)/; print $1'
echo "merging"
paste ${inFile}.chroms ${inFile}.geneids > ${inFile}.geneids_chroms.txt

echo "cleanup"
