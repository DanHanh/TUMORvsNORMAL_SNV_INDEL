
##for 001-001_001-002
# Path <- "./../BM/Annotation/104813-001-001_vs_104813-001-002/VEP/"
# InputFile <- "Mutect2_filtered_104813-001-001_vs_104813-001-002_VEP.ann.vcf.gz"
# outputFile <- "Mutations_Mutect2_filtered_104813-001-001_vs_104813-001-002_VEP.ann.xlsx"
# 
# tumorColumnNameInVcf <- 'X104813.001.001'
# normalColumnNameInVcf <- 'X104813.001.002'

##for 001-003_001-004


args = commandArgs(trailingOnly = TRUE)


project <- args[1]
projectName <- args[2]
cancerType <- args[3]
libPath <- "Rlibrary"

#Path <- paste0("./results/sarek/",project,"/Annotation/",projectName,"/VEP/")

Path <- paste0("./", project, "/Annotation/",projectName,"/VEP/")
InputFile <- list.files(Path, pattern = "Mutect2_.*\\vcf.gz$")
outputFile <- paste0("Mutations_Mutect2_",projectName, "_VEP.ann.xlsx")
cancerGeneConsortium  <- read.table('./Data/Census_allTue Dec 28 09_07_20 2021.tsv', sep = "\t", header = TRUE)
OncoKB <-read.table('./Data/cancerGeneList2.csv', sep = "\t", header = TRUE) # https://www.oncokb.org/cancerGenes
DriverDBv3 <- read.table('./Data/mutation_download_tab.txt', sep = "\t", header = TRUE) # http://driverdb.tms.cmu.edu.tw/download

## remove later
#Path <- paste0("./../results/sarek/","104202-001-003.tsv","/Annotation/","104202-001-003","/VEP/")
#InputFile <- list.files(Path, pattern = "\\vcf.gz$")
#outputFile <- paste0("Mutations_","104202-001-003", "_VEP.ann.xlsx")
#cancerGeneConsortium  <- read.table('./../Census_allTue Dec 28 09_07_20 2021.tsv', sep = "\t", header = TRUE)







library(vcfR, lib.loc = libPath)
data <- read.vcfR(paste0(Path,InputFile))


meta_info <- data@meta
vcf <- data.frame(data@fix)
vcf_gt <- data.frame(data@gt)

# get tumor sample name in vcf
library(stringr, lib.loc = libPath)
tumorColumnNameInVcf = data.frame(str_extract(meta_info, '^##tumor_sample=.*'))
tumorColumnNameInVcf <- tumorColumnNameInVcf[!is.na(tumorColumnNameInVcf)]
tumorColumnNameInVcf <- paste0("X",str_replace(tumorColumnNameInVcf , '^##tumor_sample=' , ''))
tumorColumnNameInVcf <- str_replace_all(tumorColumnNameInVcf , '-', '.')

# get normal sample name in vcf
library(stringr , lib.loc = libPath)
normalColumnNameInVcf = data.frame(str_extract(meta_info, '^##normal_sample=.*'))
normalColumnNameInVcf <- normalColumnNameInVcf[!is.na(normalColumnNameInVcf)]
normalColumnNameInVcf <- paste0("X",str_replace(normalColumnNameInVcf , '^##normal_sample=' , ''))
normalColumnNameInVcf <- str_replace_all(normalColumnNameInVcf , '-', '.')



#normalColumnNameInVcf <- ''


Df_1 <- data.frame(cbind(vcf[,1:7]))

Df_1 <- cbind(Df_1, vcf_gt)


Df_1 <-cbind( Df_1, data.frame('AS_FilterStatus' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(AS_FilterStatus=))[^;]+')))))
Df_1 <-cbind( Df_1, data.frame('AS_SB_TABLE' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(AS_SB_TABLE=))[^;]+')))))
Df_1 <-cbind( Df_1, data.frame('DP' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(DP=))[^;]+')))))
Df_1 <-cbind( Df_1, data.frame('ECNT' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(ECNT=))[^;]+')))))
Df_1 <-cbind( Df_1, data.frame('GERMQ' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(GERMQ=))[^;]+')))))
Df_1 <-cbind( Df_1, data.frame('MBQ' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(MBQ=))[^;]+')))))
Df_1 <-cbind( Df_1, data.frame('MFRL' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(MFRL=))[^;]+')))))
Df_1 <-cbind( Df_1, data.frame('MMQ' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(MMQ=))[^;]+')))))
Df_1 <-cbind( Df_1, data.frame('MPOS' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(MPOS=))[^;]+')))))
Df_1 <-cbind( Df_1, data.frame('NALOD' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(NALOD=))[^;]+')))))
Df_1 <-cbind( Df_1, data.frame('NLOD' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(NLOD=))[^;]+')))))
Df_1 <-cbind( Df_1, data.frame('POPAF' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(POPAF=))[^;]+')))))
Df_1 <-cbind( Df_1, data.frame('TLOD' = unlist(lapply(vcf$INFO, function(x) str_extract(x, '(?<=(TLOD=))[^;]+')))))





# get INFO lines from metadata informations
a <- str_extract(meta_info, '^##INFO.*')
a <- a[!is.na(a)] # removes all non INFO lines
# last ifno strings contains header of FORMAT information
InfoHeader2 <- str_extract(a[length(a)], '(?<=(Format: )).*(?=(">))') # from last INFO line get FORMAT header
header_names <- unlist(strsplit(InfoHeader2, split = "|", fixed=TRUE)) # split header 



# create data frame where all data is collected together 
library(dplyr, lib.loc = libPath)
totalNbrColumns <- (length(names(Df_1)) + length(header_names))
allColNames <- c(names(Df_1) , header_names)

finalDf <- as.data.frame(matrix("", nrow=0, ncol=totalNbrColumns))
colnames(finalDf) <- allColNames
#finalDf <- data.frame()
####################################################

# reformate List of first 22 elements, so we have a list of vectors
Df_1List <- split(Df_1, seq(nrow(Df_1)))
Df_1List <- lapply(Df_1List, function(x) as.vector(unlist(x)))

finalList <- vcf
# extract part of string that contains INFO part
finalList <- sapply(finalList$INFO,function(x) str_extract(x, '(?<=(CSQ=)).*'))
# split string if multiple entries per annotation
finalList <- sapply(finalList, function(x) strsplit(x, split=",", fixed = TRUE ))
# add white space behind last "|" so it also splits apart the last "|" in string. (otherwise different length lists if last entry behind "|" empty)
finalList <- sapply(finalList, function(x) paste0(x," "))
# Split information string by each "|"
finalList <- sapply(finalList, function(x) sapply(x, function(x) strsplit(x, split= "|", fixed = TRUE)))
# add information previously extracted in Df_1 to the list
finalList <- mapply(function(x,y) sapply(x, function(x) c(y,x) ) , finalList, Df_1List)
# change list names to 1,2,3,4 .. (names irelevant, but can get to long for data frame generation (colnames not more than 1000 bits))
names(finalList) <-  seq(1, length(finalList))
# create data frame and transpose for correct orientation
finalDf <- t(data.frame(finalList))
# name columns and rows
colnames(finalDf) <- allColNames
rownames(finalDf) <- 1:nrow(finalDf)
finalDf <- data.frame(finalDf)
###################################################

if (FALSE){ # not used anymore

# loop through each lane, extract FORMAT information, split information and add header names. Afterwards add one lane per annotation (can have multiple annotations per location)
  Length_info <- length(vcf$INFO)
  for (i in 1:length(vcf$INFO)){
    print(paste0(i, " of ", Length_info))
    string = vcf$INFO[i]
    string = str_extract(string, '(?<=(CSQ=)).*')
    string = strsplit(string, split=",", fixed = TRUE )
    string <- lapply(string,function(x) paste0(x," ")) # add empty space after each annotation. This way stringsplit adds element after final | (otherwise only if not empty!)
    string = lapply(string[[1]], function(x) strsplit(x, split= "|", fixed = TRUE))
    for (j in string){
      tmp_df <- t(j[[1]])
      colnames(tmp_df) <- header_names
      finalDf <- bind_rows(finalDf, cbind(Df_1[i,], tmp_df))
    }
  }
}
###########################################################################
## Test if in CGC
allCGC <-lapply(paste(cancerGeneConsortium$Gene.Symbol,cancerGeneConsortium$Synonyms, sep =","), function(x) unlist(strsplit(x, ",")))
allCGCmerged <- unlist(lapply(paste(cancerGeneConsortium$Gene.Symbol,cancerGeneConsortium$Synonyms, sep=","), function(x) unlist(strsplit(x, ","))))

InCGCindex <- which(finalDf$SYMBOL %in% allCGCmerged) # Variances in CGC, Indexes

finalDf$CGC_gene <- ""
finalDf$CGC_location <- ""

for(i in InCGCindex){
  gene_name <- finalDf$SYMBOL[i]
  TF <- c()
  for (cgc in allCGC){
    TF <- c(TF, any(cgc %in% gene_name))
  }
  finalDf$CGC_gene[i] <- cgcGenes <- paste0(cancerGeneConsortium$Gene.Symbol[which(TF)], collapse= ",")
  finalDf$CGC_location[i] <- paste0(cancerGeneConsortium$Genome.Location[which(TF)], collapse= ",")
}


# Test if in OncoKB

allOncoKB <-lapply(paste(OncoKB$Hugo.Symbol,OncoKB$Gene.Aliases, sep =","), function(x) trimws(unlist(strsplit(x, ","))))
allOncoKBmerged <- unlist(lapply(paste(OncoKB$Hugo.Symbol,OncoKB$Gene.Aliases, sep=","), function(x) trimws(unlist(strsplit(x, ",")))))

InOncoKBindex <- which(finalDf$SYMBOL %in% allOncoKBmerged) # Variances in CGC, Indexes

finalDf$OncoKB_gene <- ""

for(i in InOncoKBindex){
  gene_name <- finalDf$SYMBOL[i]
  TF <- c()
  for (onco in allOncoKB){
    TF <- c(TF, any(onco %in% gene_name))
  }
  finalDf$OncoKB_gene[i] <- paste0(OncoKB$Hugo.Symbol[which(TF)], collapse= ",")
}


# Test if in DriverDBv3 (remove white space in original File!!! BLCA e-driver -> in driver genes after SIN3B double ", ,")
DriverDBv3sub <- DriverDBv3[DriverDBv3$cancer_type_abbr == cancerType,]
allDriverDBv3sub <-lapply(DriverDBv3sub$driver_gene, function(x) trimws(unlist(strsplit(x, ","))))
allDriverDBv3submerged <- unlist(lapply(DriverDBv3sub$driver_gene , function(x) trimws(unlist(strsplit(x, ",")))))

InDriverDBv3subindex <- which(finalDf$SYMBOL %in% allDriverDBv3submerged) # Variances in DriverDBv3, Indexes

finalDf$DriverDBv3_tool <- ""

for(i in InDriverDBv3subindex){
  gene_name <- finalDf$SYMBOL[i]
  TF <- c()
  for (driv in allDriverDBv3sub){
    TF <- c(TF, any(driv %in% gene_name))
  }
  finalDf$DriverDBv3_tool[i] <- paste0(DriverDBv3sub$tool[which(TF)], collapse= ",")
}











tumor <- finalDf[,tumorColumnNameInVcf]
normal <- finalDf[,normalColumnNameInVcf]
format <- finalDf[,'FORMAT']
DEPTHcalc <- data.frame(format = format , normal = normal ,tumor = tumor)

TumorVsNormal <- data.frame(matrix(ncol = 4, nrow = length(DEPTHcalc$format)))
colnames(TumorVsNormal) <- c('Tumor Variant Allele Fraction','Normal Variant Allele Fraction', 'Tumor Depth', 'Normal Depth')

f <- lapply(format, function(x) unlist(strsplit(x, split = ":")))
n <- lapply(normal, function(x) unlist(strsplit(x, split = ":")))
t <- lapply(tumor, function(x) unlist(strsplit(x, split = ":")))


AlelleFrequencyNormal <- mapply(function(x,y) unlist(y)[match("AF",unlist(x))] ,f,n)
DepthNormal <- mapply(function(x,y) unlist(y)[match("DP",unlist(x))] ,f,n)
AlelleFrequencyTumor <- mapply(function(x,y) unlist(y)[match("AF",unlist(x))] ,f,t)
DepthTumor <- mapply(function(x,y) unlist(y)[match("DP",unlist(x))] ,f,t)

TumorVsNormal$`Normal Variant Allele Fraction` <- AlelleFrequencyNormal
TumorVsNormal$`Normal Depth` <- DepthNormal
TumorVsNormal$`Tumor Variant Allele Fraction` <- AlelleFrequencyTumor
TumorVsNormal$`Tumor Depth` <- DepthTumor

finalDf <- cbind(finalDf, TumorVsNormal)









selectColumns <- c('SYMBOL', 'Gene', 'Consequence', 'IMPACT', 'BIOTYPE', 'Feature_type','Feature','Tumor Variant Allele Fraction', 'Normal Variant Allele Fraction','Tumor Depth','Normal Depth', 'DriverDBv3_tool', 'OncoKB_gene','CGC_gene', 'CGC_location', 'VARIANT_CLASS','CHROM', 'POS', 'REF', 'ALT','FILTER', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons','EXON', 'INTRON' )
#selectColumns <- c('SYMBOL', 'Gene', 'Consequence', 'IMPACT', 'BIOTYPE', 'Feature_type','Feature','Tumor Variant Allele Fraction','Tumor Depth', 'CGC_gene', 'CGC_location', 'VARIANT_CLASS','CHROM', 'POS', 'REF', 'ALT','FILTER', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons','EXON', 'INTRON' )

selectedColumns <- finalDf[selectColumns]


selectedColumnsPASS <- selectedColumns[selectedColumns$FILTER == "PASS",]
selectedColumnsPASSinDB <- selectedColumnsPASS[((selectedColumnsPASS$CGC_location != "") | (selectedColumnsPASS$OncoKB_gene != "") | (selectedColumnsPASS$DriverDBv3_tool != "")),]
#selectedColumnsPASSinGCG <- selectedColumnsPASS[(selectedColumnsPASS$CGC_location != ""),]
#selectedColumnsPASSinOncoKB <- selectedColumnsPASS[(selectedColumnsPASS$OncoKB_gene != ""),]
#selectedColumnsPASSinDriverDBv3 <- selectedColumnsPASS[(selectedColumnsPASS$DriverDBv3_tool != ""),]

#a <- selectedColumnsPASS
#a <- a[(a$IMPACT == 'HIGH' | a$IMPACT == 'MODERATE'  ),]
#a <- unique(a$SYMBOL)
#a <- a[a != ""]
#write.table(a, "./mutect_PASS_M_H_geneList.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

################################################################################
#write xlsx files
#library(xlsx)
library(openxlsx, lib.loc = libPath)

print("save output file")

WB <- createWorkbook(outputFile)
addWorksheet(WB, "allData",gridLines = TRUE)
addWorksheet(WB, "OnlyPassedVariances",gridLines = TRUE)
addWorksheet(WB, "OnlyPassedVariancesInDB",gridLines = TRUE)

writeDataTable(WB, sheet="allData", selectedColumns)
writeDataTable(WB, sheet="OnlyPassedVariances", selectedColumnsPASS)
writeDataTable(WB, sheet="OnlyPassedVariancesInDB", selectedColumnsPASSinDB)

Headers <- createStyle(
  fontSize = 14,
  border = "TopBottom",
  borderStyle = "thick",
  valign = "center",
  textDecoration = "bold")
Bodies <- createStyle(
  border ="TopBottom")


addStyle(WB, sheet="allData", Headers, rows =1,cols=1:length(names(selectedColumns)) , gridExpand=TRUE)
addStyle(WB, sheet="OnlyPassedVariances", Headers, rows =1,cols=1:length(names(selectedColumns)) , gridExpand=TRUE)
addStyle(WB, sheet="OnlyPassedVariancesInDB", Headers, rows =1,cols=1:length(names(selectedColumns)) , gridExpand=TRUE)


addStyle(WB, sheet="allData", Bodies, rows= 2:length(selectedColumns$Gene), cols=1:length(names(selectedColumns)), gridExpand = TRUE)
addStyle(WB, sheet="OnlyPassedVariances", Bodies, rows= 2:length(selectedColumnsPASS$Gene), cols=1:length(names(selectedColumns)), gridExpand = TRUE)
addStyle(WB, sheet="OnlyPassedVariancesInDB", Bodies, rows= 2:length(selectedColumnsPASSinDB$Gene), cols=1:length(names(selectedColumns)), gridExpand = TRUE)


setColWidths(WB, sheet="allData",cols=1:length(names(selectedColumns)), widths = "auto")
setColWidths(WB, sheet="OnlyPassedVariances",cols=1:length(names(selectedColumnsPASS)), widths = "auto")
setColWidths(WB, sheet="OnlyPassedVariancesInDB",cols=1:length(names(selectedColumnsPASSinDB)), widths = "auto")

setRowHeights(WB, sheet="allData", rows= 1:length(selectedColumns$Gene), heights = 20)
setRowHeights(WB, sheet="OnlyPassedVariances", rows= 1:length(selectedColumnsPASS$Gene), heights = 20)
setRowHeights(WB, sheet="OnlyPassedVariancesInDB", rows= 1:length(selectedColumnsPASSinDB$Gene), heights = 20)


##############################
#write metadata, all columns DF and CancerGeneCensus

addWorksheet(WB, "CancerGeneCensus",gridLines = TRUE)
#addWorksheet(WB, "MetadataAllColumns",gridLines = TRUE)
#addWorksheet(WB, "AllColumnsData",gridLines = TRUE)

writeDataTable(WB, sheet="CancerGeneCensus", cancerGeneConsortium)
#writeData(WB, sheet="MetadataAllColumns", data.frame(a))
#writeDataTable(WB, sheet="AllColumnsData", finalDf)

addStyle(WB, sheet="CancerGeneCensus", Headers, rows =1,cols=1:length(names(cancerGeneConsortium)) , gridExpand=TRUE)
#addStyle(WB, sheet="AllColumnsData", Headers, rows =1,cols=1:length(names(finalDf)) , gridExpand=TRUE)

addStyle(WB, sheet="CancerGeneCensus", Bodies, rows= 2:length(cancerGeneConsortium$Gene.Symbol), cols=1:length(names(cancerGeneConsortium)), gridExpand = TRUE)
#addStyle(WB, sheet="AllColumnsData", Bodies, rows= 2:length(finalDf$CHROM), cols=1:length(names(finalDf)), gridExpand = TRUE)

setColWidths(WB, sheet="CancerGeneCensus",cols=1:length(names(cancerGeneConsortium)), widths = "auto")
#setColWidths(WB, sheet="AllColumnsData",cols=1:length(names(finalDf)), widths = "auto")

setRowHeights(WB, sheet="CancerGeneCensus", rows= 1:length(cancerGeneConsortium$Gene.Symbol), heights = 20)
#setRowHeights(WB, sheet="AllColumnsData", rows= 1:length(finalDf$CHROM), heights = 20)

saveWorkbook(WB, paste0(Path, outputFile), overwrite = TRUE)

#copy file to excel location
#file.copy(paste0(Path, outputFile), paste0("./results/excelFiles/", outputFile), overwrite=TRUE)
