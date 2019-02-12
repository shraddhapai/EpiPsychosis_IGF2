#!/bin/bash

# cat fasta files for bsmap 
inDir=/home/shraddhapai/genome_annotation/grch38/fasta/chroms
lambda=/home/shraddhapai/genome_annotation/lambda.fa

outFile=${inDir}/hg38NoAlt_Lambda.fa

declare -a chroms=({1..22} X Y M);

cat /dev/null > $outFile
for i in {0..24}; do
	echo ${chroms[i]}
	cat ${inDir}/chr${chroms[i]}.fa >> $outFile
done
cat $lambda >> $outFile

echo $cmd
`$cmd`

