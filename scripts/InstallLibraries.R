
args = commandArgs(trailingOnly = TRUE)
print(args[1])
libPath=args[1]

.libPaths(libPath)

repo <- "https://cloud.r-project.org/"

install.packages("BiocManager", lib = libPath, repos = repo)

install.packages("vcfR",lib = libPath, repos = repo,dependencies = TRUE, force = TRUE)
install.packages("stringr",lib = libPath, repos = repo,dependencies = TRUE, force = TRUE)
install.packages("dplyr",lib = libPath, repos = repo,dependencies = TRUE, force = TRUE)
install.packages("openxlsx",lib = libPath, repos = repo, dependencies = TRUE, force = TRUE)

