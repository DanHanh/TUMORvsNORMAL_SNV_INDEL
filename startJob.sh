#!/bin/bash

InputFile="$1"
outDir="$2"
cancerType="$3"
d=$(date)
mkdir RunConditions

echo "---------------------------------------------------" >> ./RunConditions/runLogs.txt
echo $d >> runLogs.txt
echo "InputFile: $1" >> ./RunConditions/runLogs.txt
echo "outDir: $2" >> ./RunConditions/runLogs.txt
echo "cancerType: $3" >> ./RunConditions/runLogs.txt
echo "---------------------------------------------------" >> ./RunConditions/runLogs.txt

sbatch --job-name="$outDir" ./scripts/job.sh "$InputFile" "$outDir" "$cancerType"
