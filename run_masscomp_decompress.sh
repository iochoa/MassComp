#!/bin/bash

# Input argument is the folder containing the masscomp compressed files that need to be decompressed

# Please modify the outputPath
# Please modify the location of the masscomp executable

outputPath='/srv/data/ruochen/bioinformatics/data/masscomp_results/decompressed/'

for inputFile in $1/*.mc; do
	name=$(basename "$inputFile")
	echo $inputFile
	#echo $name
	../../MassComp-master/src/masscomp -d $inputFile $outputPath/$name.mc.mzXML

done


