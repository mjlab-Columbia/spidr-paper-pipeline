#!/bin/bash
file=splitclusters.in #directorywheresplitclusterfilesare
SUBMISSION_NAME=qsub.sh
bamfile="/central/scratch/mblanco/pipelinePB0627/workup/alignments/fibroblast_SPIDR_5p_RNA.merged.RPM.bam"
while read line           
do           
    	echo "$line"	
	filename=$(basename $line)
	echo $filename
	newdir=$filename"splitbamminoligo1"
	directory=$(dirname $line)
	cd $directory
	mkdir $newdir
	cd $newdir

# Quality trimming with Trimmomatic
    echo \#!/bin/bash > $SUBMISSION_NAME
    echo \#SBATCH --time=16:00:00 >> $SUBMISSION_NAME  # walltime
    echo \#SBATCH --ntasks=1 >> $SUBMISSION_NAME   # number of processor cores (i.e. tasks)
    echo \#SBATCH --nodes=1 >> $SUBMISSION_NAME  # number of nodes
    echo \#SBATCH --mem-per-cpu=16G >> $SUBMISSION_NAME   # memory per CPU core
    echo \#SBATCH --cpus-per-task=8 >> $SUBMISSION_NAME #cpus per taks
    echo \#SBATCH -J "trimmomatic" >> $SUBMISSION_NAME    # job name

    #Trimmomatic Paired End
    echo "python /groups/guttman/projects/spidr/20220916/mtorqcnewconfig/scripts/threshold_tag_and_split.py -i $bamfile -o $filename".bam" -c $line --num_tags 8 -d $newdir --min_oligos 1 --proportion 0.8 --max_size 100" >> $SUBMISSION_NAME

    # submit
    sbatch $SUBMISSION_NAME 

done < $file
