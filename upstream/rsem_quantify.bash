#!/bin/bash	
	
homedir="/rds/general/user/mmoradim/home/RNAseq_SCVD19"
resultsDir="/rds/general/user/mmoradim/projects/cunnington-covid-transcriptomics/ephemeral/Users/Mahdi/RNAseq_SCVD19/rsem_quantify_output"
	
scriptsPath="$homedir/scripts/rsem_quantify"
logDir="$scriptsPath/logs"
mkdir -p $logDir	
	
inPath="/rds/general/user/mmoradim/projects/cunnington-covid-transcriptomics/ephemeral/Users/Mahdi/RNAseq_SCVD19/star_map_output"
inPathBit="Aligned.toTranscriptome.out.bam"	

names=($(ls -d $inPath/*))
	
for dir in ${names[@]};do
        shortname=`basename $dir`
       	echo $shortname
		inFile=( $(ls $dir/*$inPathBit) )
		outPathBit=`basename $inFile | sed 's/_Aligned.toTranscriptome.out.bam//'`	
        echo $outPathBit	
        outFilePrefix=$resultsDir/$outPathBit	
                	
        rsem_quantify_version="/rds/general/user/mmoradim/home/tools/RSEM-master/rsem-calculate-expression"	
                	
        reference="/rds/general/user/mmoradim/projects/cunnington-covid-transcriptomics/ephemeral/Users/Mahdi/RNAseq_SCVD19/rsem_index/hg38"
                	
        rsem_quantify_line="$rsem_quantify_version \
					--paired-end \
					--bam \
					--no-bam-output \
					--seed 12345 \
					-p 6 \
					--forward-prob 0 \
					$inFile \
					$reference \
					$outFilePrefix ;"
	
	rsem_quantify_command="rsem_quantify_$outPathBit"
	echo $rsem_quantify_line >> $rsem_quantify_command
	chmod 775 $rsem_quantify_command
        # qsub -l select=1:mem=10gb:ncpus=10 -e $logDir -N rsem_quantify\_$outPathBit -l walltime=72:00:00 $rsem_quantify_command	
      	
qsub <<< "	
#PBS -l select=1:mem=20gb:ncpus=6	
#PBS -e $logDir	
#PBS -N $shortname
#PBS -l walltime=72:00:00	
	
module load perl/5.24.0	
$scriptsPath/$rsem_quantify_command"	
	  
done;	
