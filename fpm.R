# counts <- read.delim('output/counts/featurecounts/hisat2/celegans/ensembl/release96/countsMinoverlap1strandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600DiscardMisalignedPairsNoMultimappers.txt')
readsMappedPaired <- matrix(0, nrow = 18, ncol = 4)
readsMappedPaired <- as.data.frame(readsMappedPaired)
colnames(readsMappedPaired) <- c('stage', 'bact', 'repl', 'nbReads')

stages <- c('adults', 'embryos')
bacteria <- c('hb101', 'BIGb446', 'BIGb468')
replicates <- c('repl1', 'repl2', 'repl3')

i <- 1

for (stage in stages) {
    for (bact in bacteria) {
        for (repl in replicates) {
            filePath <- paste0('output/alignment/target/celegans/ref/ensembl/release96/query/', stage, '/', bact, '/', repl, '/hisat2/alnSortedByCoordStatsReadsMappedAndPaired.txt')
            readsMappedPaired[i, ] <- c(stage, bact, repl, read.table(filePath,  header = FALSE)[1,1])
            i <- i + 1
        }
    }
}

write.csv(readsMappedPaired, 'output/alignment/target/celegans/ref/ensembl/release96/query/readsMappedPaired.csv', row.names = FALSE, quote = FALSE)
