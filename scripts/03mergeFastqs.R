### Merge fastq files together that come from the 
### same sample and have reads in the same direction
### but have been produced in different lanes 
### (The sequencing technician from the Gurdon Institute showed me the 2 lanes on the machine)

library(stringr)

### This script takes 4 arguments
# stage : C. elegans synchronised population stage
# bact : bacteria the worms were exposed to
# repl : replicate
# readOrientation : read orientation (forward or reverse orientation)

stage <- snakemake@wildcards$stage
bact <- snakemake@wildcards$bact
repl <- snakemake@wildcards$repl
readOrientation <- snakemake@wildcards$readOrientation
fastqFiles <- dir(paste0('output/reads/rna/illumina/celegans/raw/', stage, '/', bact, '/', repl, '/'))
fastqFilesPaths <- paste0('output/reads/rna/illumina/celegans/raw/', stage, '/', bact,
                   '/', repl, '/', 
                   fastqFiles[stringr::str_detect(fastqFiles, readOrientation)])


system(paste0('cat ', fastqFilesPaths[1], ' ', fastqFilesPaths[2],
       ' > output/reads/rna/illumina/celegans/merged/', stage,
       '/', bact, '/', repl, '/R', readOrientation, '.fq.gz'))
