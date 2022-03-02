#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --mem=128000
#SBATCH --job-name=BLCa128
#SBATCH --cpus-per-task=16

module load Workspace_Home;
module load Nextflow/20.10.0

InputFile="$1"
outDir="$2"
cancerType="$3"


nextflow run nf-core/sarek -r 2.7.1 --tools Manta,MSIsensor,Strelka,Mutect2,VEP,snpEff --save_bam_mapped true --outdir "$outDir" --input "$InputFile" --genome GRCh38 -profile singularity --step mapping --trim_fastq true

mkdir ./${outDir}/excelFiles


Path=${outDir}'/Annotation/'$( ls ./${outDir}/Annotation/ | grep '.*_vs_.*' )'/VEP/'

projectName=$( ls ./${outDir}/Annotation/ | grep '.*_vs_.*' )

Rscript ./scripts/Mutect2ProcessAnnVcfFiles2.0.R "$outDir" "$projectName" "$cancerType"
Rscript ./scripts/StrelkaSNVsProcessAnnVcfFiles2.0.R "$outDir" "$projectName" "$cancerType"
Rscript ./scripts/StrelkaIndelProcessAnnVcfFiles2.0.R "$outDir" "$projectName" "$cancerType"

#copy xlsx files to excelfile location
cp ${Path}/*.xlsx ./${outDir}/exelFiles/
