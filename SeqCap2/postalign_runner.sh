#!/bin/bash

# wrapper to postalign.sh whose sole purpose is to launch log files
# for each sample
declare -a SAMPLE=(10Redoneg 10Redopos 20Redoneg 20Redopos 23Redoneg 23Redopos 24Redoneg 24Redopos 25Redoneg 25Redopos 27neg 27pos 2Redoneg 2Redopos 33neg 33pos 37neg 37pos 4neg 4pos 58neg 58pos 66Redoneg 66Redopos 69neg 69pos 70neg 70pos 74neg 74pos 75neg 75pos 79neg 79pos 80neg 80pos 81neg 81pos 87neg)

alignDir=/scratch/g/gbader/spai/Epigenetics/NARSAD/input_files/SeqCapEpi2/LABV_20180618_Capture/align

jobdir="postalignjobs"
mkdir -p $jobdir

arraylen=${#SAMPLE[@]}
echo $arraylen
curd=`pwd`
for i in {1..2}; do
	sampName=${SAMPLE[$i]}
	echo "**$sampName**"
	logFile=${alignDir}/${sampName}.postalign
	jobfile=${jobdir}/${sampName}.postalign.sh
	echo $jobfile
	cat /dev/null > $jobfile
	echo "#!/bin/bash" >> $jobfile
	echo "#SBATCH --nodes=1" >> $jobfile
	echo "#SBATCH --cpus-per-task=30" >> $jobfile
	echo "#SBATCH --time=05:00:00" >> $jobfile
	echo "#SBATCH --job-name ${sampName}_postalign" >> $jobfile
	echo "#SBATCH --output ${alignDir}/${sampName}_align.txt" >> $jobfile
	cmd="${curd}/postalign.sh $sampName $alignDir" 
	echo $cmd >> $jobfile
	chmod u+x $jobfile
	sbatch $jobfile
done

