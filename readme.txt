# to run necessary R packages run(not with sbatch):
./installDependencies.sh

# to run tumor vs normal run:
# cancer Types (for Driver3DBv3 ): {"ACC"  "BLCA" "BRCA" "CESC" "CHOL" "COAD" "DLBC" "ESCA" "GBM" #"HNSC" "KICH" "KIRC" "KIRP" "LAML" "LGG"  "LIHC" "LUAD" "LUSC" "MESO" "OV"   "PAAD" "PCPG" "PRAD" #"READ" "SARC" "SKCM" "STAD" "TGCT" "THCA" "THYM" "UCEC" "UCS"  "UVM"}

./startJob.sh <exampleInputFile.tsv> <outDir> <cancerType>
