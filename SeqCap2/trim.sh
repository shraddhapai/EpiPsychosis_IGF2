#!/bin/bash
# trimmomatic to remove low-quality reads from the fastq files.

# trimmomatic user manual and helpful tutorials:
# http://www.usadellab.org/cms/?page=trimmomatic
# http://rnaseq.agbioinfo.utk.edu/index.php/Trimming_adapters_and_low_quality_reads_(Trimmomatic)

TRIMMO=/home/shraddhapai/software/Trimmomatic-0.36/trimmomatic-0.36.jar
trimAdapter=/home/shraddhapai/software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa
FASTQC=/home/shraddhapai/software/FastQC/fastqc

rootDir=/home/shraddhapai/Epigenetics/NARSAD/input_files/SeqCapEpi2/LABV_20180618_Capture
logdir=${outdir}

indir=${rootDir}/fastq
outdir=${rootDir}/trimmed
unpairdir=${outdir}/unpaired
mkdir -p $outdir
mkdir -p $unpairdir

echo $indir
for r1 in ${indir}/*R1*gz;do
	baseF=`basename $r1`
	echo $baseF
	r2=${baseF/R1/R2}
	echo $r2

	out1=${r1/.fastq.gz/.trim.fastq.gz}
	out1=`basename $out1`
	out2=${r2/.fastq.gz/.trim.fastq.gz}

	up1=${out1/trim/unpaired}
	up2=${out2/trim/unpaired}

	java -jar $TRIMMO PE -threads 6 -phred33 \
		$r1 ${indir}/${r2} \
		${outdir}/${out1} ${outdir}/${up1} \
		${outdir}/${out2} ${outdir}/${up2} \
		ILLUMINACLIP:${trimAdapter}:2:30:10 \
		LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:75  &>  ${outdir}/${baseF}.log
done

# count trimmed file reads
outfile=trimmed_fastq_countLines.log
cat /dev/null > $outfile
for f in ${outdir}/*trim.fastq.gz;do
	baseF=`basename $f`
	echo $baseF
	nl=`zcat $f | wc -l | awk '{rd=$1 / 4; print rd}'`
	echo -e "$baseF\t$nl" >> $outfile
	###$FASTQC $f --outdir=$outdir --extract
done
