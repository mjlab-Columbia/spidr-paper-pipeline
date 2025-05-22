#!/bin/bash

proteinfile="proteins.in"

#merge bam files
while read line
do
	rm mergesam.in
	echo $line
	#echo "*$line".bam"" 
	find . -wholename "*minoligo3*$line".bam"" -type f > mergesam.in
	samtools merge -f -b mergesam.in $line".merged.minoligo3.bam"

done < $proteinfile
