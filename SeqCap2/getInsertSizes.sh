#!/bin/bash

bamDir=/home/shraddhapai/Epigenetics/NARSAD/output_files/SeqCap2/clean_bam
samtools=~/software/samtools-1.5/samtools

cat /dev/null > insert_stats.txt
for fName in ${bamDir}/*CLEAN.bam;do
	baseF=`basename $fName`
	echo $baseF
	val=`$samtools view $fName |  sort -u | awk 'BEGIN {tot=0; ctr=0; ss=0} 
		{ tmp=length($10);  tot+=tmp; ctr+=1; ss+=tmp^2;} 
	 END { avg=tot/ctr; 
		sd=sqrt((ss/NR)-(tot/NR)^2);
		print avg"\t "sd; 
	} '`
	echo -e "$baseF\t$val" >> insert_stats.txt
	#mu=`grep "insert size average" tmp.txt | cut -f 3`
	#sd=`grep "insert size standard" tmp.txt | cut -f 3`
	#echo -e "$baseF\t$mu\t$sd" >> insert_stats.txt
done
