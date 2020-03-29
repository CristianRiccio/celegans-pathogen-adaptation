library(tidyr)

# Read in the featurecounts file
counts <- read.table('output/HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt',
                     header = TRUE)

