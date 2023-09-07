#!/bin/bash

homedir="/rds/general/user/mmoradim/home/RNAseq_SCVD19"
genomeDir="/rds/general/user/mmoradim/projects/cunnington-covid-transcriptomics/ephemeral/Users/Mahdi/RNAseq_SCVD19/star_index"
resultsDir="/rds/general/user/mmoradim/projects/cunnington-covid-transcriptomics/ephemeral/Users/Mahdi/RNAseq_SCVD19/star_map_output"

scriptsPath="$homedir/scripts/star_map"
logDir=$scriptsPath"/logs"
mkdir -p $logDir

inPath="/rds/general/user/mmoradim/projects/cunnington-covid-transcriptomics/ephemeral/Users/Mahdi/RNAseq_SCVD19/trimmed_data"
inPathBit="_*.trimmed.fastq.gz"

names=($(ls -d $inPath/*))

for dir in ${names[@]};do
	shortname=`basename $dir`
	echo $shortname
	mkdir $resultsDir/$shortname
	files=( $(ls $dir/*$inPathBit) )
	inFile1=${files[0]}
	inFile2=`echo $inFile1 | sed 's/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz/'`

        star_map_line="STAR \
				--genomeDir $genomeDir \
                --readFilesIn $inFile1 $inFile2 \
                --readFilesCommand zcat \
                --outFileNamePrefix $resultsDir/$shortname/$shortname"_" \
                --outFilterType BySJout \
                --outSAMtype BAM SortedByCoordinate \
                --outWigType wiggle \
				--outWigNorm None \
                --quantMode TranscriptomeSAM ;"
                
	star_map_command="star_map_$shortname"
	echo $star_map_line >> $star_map_command
	chmod 775 $star_map_command
	# qsub -l select=1:mem=40g:ncpus=2 -e $logDir -N $shortname -l walltime=72:00:00 $star_map_command

qsub <<< "
#PBS -l select=1:mem=40gb:ncpus=2
#PBS -e $logDir
#PBS -N $shortname
#PBS -l walltime=72:00:00

module load star/2.7.1a
$scriptsPath/$star_map_command"

done;