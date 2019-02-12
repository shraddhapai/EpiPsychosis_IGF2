#!/bin/bash

plink=/home/shraddhapai/software/plink/plink
snpFile=/home/shraddhapai/Epigenetics/NARSAD/anno/EPIC_array_SNPs.txt
genoF=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/SUS19399.hg19.sorted_-clean

outDir=`dirname $genoF`
$plink --bfile $genoF --extract $snpFile --recode A --out ${genoF}.EPICsnps
