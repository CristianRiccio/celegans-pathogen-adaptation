### Download the data from here
### https://basespace.illumina.com/run/124171047/HS678
### to ~/Documents/phd/nick/input/rawData/


### No need to remove indices, already done by the demultiplexing

### No need to remove adapter, the read should not reach it

### N2 fed HB101, BIGb446 and BIGb468

### Create a conda environment for Nick's analysis
# conda create --yes --name nick bioconductor-all bioconductor-annotationdbi bioconductor-biomart bioconductor-deseq2 bioconductor-topgo bioconductor-tximport cutadapt imagemagick kallisto r-pheatmap r-stringr r-yaml samtools snakemake star subread
# conda create --yes --name nick2 python=2.7

# withr::with_libpaths(new = 'output/software/r/package/installation/',
#                      devtools::install_github('pachterlab/sleuth', 'hadley/colformat'))

library(sleuth, lib.loc = 'output/software/r/package/installation/')
# source("https://bioconductor.org/biocLite.R")
# biocLite("ALL")

library(DESeq2)
library(stringr)
library(sleuth)

# ### Rename the folder names
# oldNames <- dir('input/rawData')
# newNames <- stringr::str_replace(oldNames, '-ds[.].*', replacement = '')
# for (i in 1:length(oldNames)) {
#   system(paste0('mv input/rawData/', oldNames[i],
#                 ' input/rawData/', newNames[i]))
# }

# ### Uncompress all the FASTQ files only if necessary
# folderNames <- dir('input/rawData/')
# for (i in 1:length(folderNames)) {
#   fileNames <- dir(paste0('input/rawData/', folderNames[i]))
#   ### Select the compressed file names
#   compressed <- fileNames[substr(fileNames, nchar(fileNames) - 2, nchar(fileNames)) == '.gz']
#   ### If there are no compressed files, jump to the next iteration
#   if (length(compressed) == 1) next
#   
#   for (j in 1:length(compressed)) {
#     system(paste0('gunzip input/rawData/', 
#                   folderNames[i], '/', compressed[j]))
#   }
# }

