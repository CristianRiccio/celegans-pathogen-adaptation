# Calculate FPKMs (
# =fragment per kilobase of CDS per million mapped and paired fragments)

# Assign snakemake wildcards to R variables
# HSrunNumber <- 678
HSrunNumber <- snakemake@wildcards$HSrunNumber
# database <- 'ensembl'
database <- snakemake@wildcards$database
# databaseRelease <- 96
databaseRelease <- snakemake@wildcards$databaseRelease
# strandedness <- 'unstranded'
strandedness <- snakemake@wildcards$strandedness

# Load R packages
library(dplyr)
library(stringr)

# Read in the CSV file containing the number of mapped and 
# paired reads for each sample
# path2readsMappedPaired <- paste0('output/HS', HSrunNumber, 
#                                  '/alignment/target/celegans/ref/', database, 
#                                  '/release', databaseRelease,
#                                  '/query/readsMappedPaired.csv')
path2readsMappedPaired <- snakemake@input$readsMappedPaired
readsMappedPaired <- read.csv(path2readsMappedPaired)

# Read in the featurecounts file
# path2counts <- paste0('output/HS', HSrunNumber, 
#                       '/counts/featurecounts/hisat2/celegans/',
#                        database, '/release', databaseRelease, 
#                       '/minoverlap1', strandedness, 'FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsRenamed.txt')
path2counts <- snakemake@input$counts
counts <- read.table(path2counts, header = TRUE)

# Create a data frame that will hold the fpms
fpm <- counts

for (i in 7:ncol(fpm)) {
  samplei <- stringr::str_split(colnames(counts)[i], '[.]')[[1]][c(10, 12, 13)]
  wormStraini <- samplei[1]
  bacteriumi <- samplei[2]
  replicatei <- str_split(samplei[3], '')[[1]][5]  %>% as.numeric()
  
  nbReadsi <- dplyr::filter(readsMappedPaired, wormStrain == wormStraini, bacterium == bacteriumi, replicate == replicatei) %>% 
    dplyr::select(nbReads) %>% as.numeric()
  fpm[, i] <- counts[, i]/((nbReadsi/2)/1e6) # number of reads / 2 to get number of fragments
}

# Transform the FPMs into FPKMs
fpkm <- fpm

for (i in 7:ncol(fpm)) {
  fpkm[, i] <- fpm[, i]/(fpm$Length/1000)
}

# Check manually for one gene in one sample that FPKM is correct
# Gene WBGene00015153 has length 1327
# sample output.HS755.alignment.target.celegans.ref.ensembl.release96.query.n2.embryos.hb101.repl1.hisat2.alnSortedByCoord.bam
# counts[3, 7] has 7 counts
# N2 embryos HB101 replicate 1 has 46243858 mapped and paired reads
# fpkm[3, 7] == (counts[3, 7]/((46243858/2)/1e6))/(1327/1000) # TRUE

# path2fpkms <- paste0('output/HS', HSrunNumber, 
#                      '/counts/featurecounts/hisat2/celegans/',
#                      database, '/release', databaseRelease, 
#                      '/minoverlap1', strandedness, 
#                      'FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersFPKMs.txt')
path2fpkms <- snakemake@output$fpkm
print(path2fpkms)
write.csv(fpkm, path2fpkms, row.names = FALSE, quote = FALSE)