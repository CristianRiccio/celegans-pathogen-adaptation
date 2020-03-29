# Create a CSV table containing the number of mapped and paired reads
# for each sample

# Load the dplyr package
library(dplyr)

# Import the samples information
# samples <- read.delim('samples/samples.tsv')
samples <- read.delim(snakemake@input$samples)

# Assign the wildcards to variables
# HSrunNumberi <- '755'
HSrunNumberi <- snakemake@wildcards$HSrunNumber
# databasei <- 'ensembl'
databasei <- snakemake@wildcards$database
# databaseReleasei <- '96'
databaseReleasei <- snakemake@wildcards$databaseRelease

# Select the samples for the HS number of interest
samplesRun <- dplyr::filter(samples, HSrunNumber == HSrunNumberi)

# Initialise a column that will contain the number of mapped and paired reads
samplesRun$nbReads <- -1

# Loop through each sample, retrieve the number of mapped 
# and paired reads and add it to the table
for (i in 1:nrow(samplesRun)) {
  
  wormStraini <- samplesRun$wormStrain[i]
  stagei <- samplesRun$stage[i]
  bacteriumi <- samplesRun$bacterium[i]
  replicatei <- samplesRun$replicate[i]
  
  path2mappedAndPairedi <- paste0('output/HS', HSrunNumberi, 
                                  '/alignment/target/celegans/ref/', 
                                  databasei, '/release', databaseReleasei, 
                                  '/query/', wormStraini, '/', stagei, 
                                  '/', bacteriumi, '/repl', replicatei, 
                                  '/hisat2/alnSortedByCoordStatsReadsMappedAndPaired.txt')
  mappedAndPairedi <- read.table(path2mappedAndPairedi)[1, 1]
  samplesRun$nbReads[i] <- mappedAndPairedi
}


# Assign the path to the CSV table to a variable
# csvTable <- paste0('output/HS', HSrunNumberi, 
#                    '/alignment/target/celegans/ref/',
#                    databasei, '/release', databaseReleasei,
#                    '/query/readsMappedPaired.csv')
csvTable <- snakemake@output$csvTable

# Write out the dataframe containing the samples info + 
# number of mapped and paired reads to a CSV file
write.csv(samplesRun, csvTable, row.names = FALSE, quote = FALSE)
  
  