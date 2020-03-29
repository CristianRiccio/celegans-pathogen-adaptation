# Subset the counts table so that they only contain counts and no metadata (except gene ID as row names)

# countdata <- read.table('output/HS678HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt', 
#                         header = TRUE, row.names = 1)
countdata <- read.table(snakemake@input$counts, 
                        header = TRUE, row.names = 1)

# Remove first five columns (chr, start, end, strand, length)
countdata2 <- countdata[, 6:ncol(countdata)]

# write.csv(countdata2, 'output/HS678HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsOnly.csv',
#           quote = FALSE)
write.csv(countdata2, snakemake@output$countsOnly,
          quote = FALSE)