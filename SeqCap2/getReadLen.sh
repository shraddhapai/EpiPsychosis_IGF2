#!/bin/bash

# look at first read in fastq file to get read length

fqDir= # insert here
# 
zcat 76-2re_S8_L001_R1_001.trim.fastq.gz | head -n 8 | awk '{if(NR%4==2) print length($1)}' > input.readslength.txt
