# Plot the number of mapped and paired reads for each sample

# Load the 'ggplot2' package
library(ggplot2)

# Assign the wildcards to variables
HSrunNumber <- snakemake@wildcards$HSrunNumber
database <- snakemake@wildcards$database
databaseRelease <- snakemake@wildcards$databaseRelease

# Read in the CSV file containing the number of mapped and 
# paired reads for each sample
# readsMappedPaired <- read.csv('output/HS755/alignment/target/celegans/ref/ensembl/release96/query/readsMappedPaired.csv')
readsMappedPaired <- read.csv(snakemake@input$readsMappedPaired)

# Transform the replicate column to a factor
readsMappedPaired$replicate <- as.factor(readsMappedPaired$replicate)

# Plot a barplot of the numbers of mapped and paired reads
ggplot2::ggplot(readsMappedPaired, ggplot2::aes(x = replicate, y = nbReads, fill = replicate)) +
  ggplot2::geom_bar(position = 'dodge', stat = 'identity') +
  ggplot2::ggtitle('Number of mapped and paired reads') +
  ggplot2::scale_y_continuous(labels = function(x) x/1e6) +
  ggplot2::ylab('Number of reads / million') +
  ggplot2::facet_grid(~ bacterium + wormStrain)
# ggplot2::ggsave('output/HS755/images/bar/celegans/ref/ensembl/release96/nbMappedPairedReads.pdf', width = 7, height = 7)
ggplot2::ggsave(snakemake@output$plot, width = 7, height = 7)

