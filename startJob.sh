#!/bin/bash


InputFile="$1"
outDir="$2"
cancerType="$3"
sbatch ./scripts/job.sh "$InputFile" "$outDir" "$cancerType"
