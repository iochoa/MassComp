#!/bin/bash

# The input argument is a folder containing the mxXML files to be compressed
# Please modify the outputPath accordingly, as well as the location of the masscomp executable

outputPath='/srv/data/ruochen/bioinformatics/data/masscomp_results/MSV000080905/Plate4/'

for inputFile in $1/*.mzXML; do
	name=$(basename "$inputFile")
	echo $inputFile
	#echo $name
	../../MassComp-master/src/masscomp -c $inputFile $outputPath/$name.mc

done


