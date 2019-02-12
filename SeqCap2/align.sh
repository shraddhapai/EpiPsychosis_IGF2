#!/bin/bash

# align BSseq reads with BSmap
# BSMap params: https://github.com/genome-vendor/bsmap/blob/master/README.txt
BSMAP=/home/g/gbader/spai/software/bsmap-2.74/bsmap
grch38_fa=/home/g/gbader/spai/genome_annotation/grch38/fasta/hg38NoAlt_Lambda.fa 
numCores=24; # Dell workstation
rootDir=/scratch/g/gbader/spai/Epigenetics/NARSAD/input_files/SeqCapEpi2/LABV_20180618_Capture

rawDir=${rootDir}/trimmed
alignDir=${rootDir}/align
mkdir -p $alignDir

jobdir="alignjobs"
mkdir -p $jobdir

for r1 in ${rawDir}/*R1*.trim.fastq.gz; do
 	baseF=`basename $r1`
    echo $baseF
    r2Base=${baseF/R1/R2}
    echo $r2Base
	r2=${rawDir}/${r2Base}
	
	#oF=`basename $r1 _R1_001.trim.fastq.gz`
	oF=`basename $r1 _R1_001.trim.fastq.gz`
	oF=${alignDir}/${oF}.bam
	echo "Outfile is $oF"
	# -r 0 : report only unique hits/pairs only
	# -s 16 : seed length
	# -n 1: Map R1 and R2 to both + and - strand (Cokus protocol)
	cmd="$BSMAP -r 0 -s 16 -n 1 -a $r1 -b $r2 -d $grch38_fa -p 30 -o $oF"

	jobfile=${jobdir}/${baseF}.align.sh
	cat /dev/null > $jobfile
	echo "#!/bin/bash" >> $jobfile
	echo "#SBATCH --nodes=1" >> $jobfile
	echo "#SBATCH --cpus-per-task=30" >> $jobfile
	echo "#SBATCH --time=00:40:00" >> $jobfile
	echo "#SBATCH --job-name ${baseF}_align" >> $jobfile
	echo "#SBATCH --output ${alignDir}/${baseF}_align.txt" >> $jobfile
	echo $cmd >> $jobfile
	chmod u+x $jobfile
	cd $jobdir
	#sbatch ${baseF}.align.sh
	cd ..
done
