#!/bin/bash

# wrapper to postalign.sh whose sole purpose is to launch log files
# for each sample
alignDir=/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi/output/align
declare -a SAMPLE=(31-1re 31-2re 31-3re 36re 41re 76-1re 76-2re 76-3re 77re 88re 90re 95re);

#alignDir=/home/shraddhapai/Epigenetics/NARSAD/SeqCapEpi/output/test

for i in 0; do
	sampName=${SAMPLE[$i]}
	logFile=${alignDir}/${sampName}.postalign
	./postalign.sh $sampName $alignDir 2>&1 | tee -a ${logFile}.log
done

