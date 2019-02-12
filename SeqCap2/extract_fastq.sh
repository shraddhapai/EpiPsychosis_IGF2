#!/bin/bash

srcFile=/media/shraddhapai/NARSAD_EpiPsychosis/SeqCap2/LABV_20180618_Capture.tar.gz
outDir=/home/shraddhapai/Epigenetics/NARSAD/input_files/sendToRoche

# 25 
tar xvfz ${srcFile} LABV_20180618_Capture/flowcell1/25Redopos_L000_R2_001.fastq.gz -C $outDir
tar xvfz ${srcFile} LABV_20180618_Capture/flowcell1/25Redopos_L000_R1_001.fastq.gz -C $outDir
# 81
tar xvfz ${srcFile} LABV_20180618_Capture/flowcell1/81pos_L000_R1_001.fastq.gz -C $outDir
tar xvfz ${srcFile} LABV_20180618_Capture/flowcell1/81pos_L000_R2_001.fastq.gz -C $outDir
# 80neg
tar xvfz ${srcFile} LABV_20180618_Capture/flowcell2/80neg_L000_R1_001.fastq.gz -C $outDir
tar xvfz ${srcFile} LABV_20180618_Capture/flowcell2/80neg_L000_R2_001.fastq.gz -C $outDir

# 95re
cp /media/shraddhapai/NARSAD_EpiPsychosis/SeqCap/SeqCapEpiChoice/Samples/95re/Files/*L001* $outDir/.

