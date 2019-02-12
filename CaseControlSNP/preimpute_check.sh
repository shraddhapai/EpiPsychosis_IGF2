#!/bin/bash

# recommended by Mich Impute Server because without it imputation fails

PLINK=/home/shraddhapai/software/plink/plink

inFile=/home/shraddhapai/Epigenetics/NARSAD/output_files/CaseControlSNP/plinkQC/hg19/SUS19399.hg19
# http://www.well.ox.ac.uk/~wrayner/tools/#Checking
CHECKER=/home/shraddhapai/software/preimpute_check/HRC-1000G-check-bim-v4.2.pl
KGP_REF=/home/shraddhapai/software/preimpute_check/1000GP_Phase3_combined.legend

echo "* write frequency file"
head ${inFile}.bim
#$PLINK --bfile $inFile --freq 
echo "* running checker"
#perl $CHECKER -b ${inFile}.bim -f plink.frq -r $KGP_REF -g -p ALL


