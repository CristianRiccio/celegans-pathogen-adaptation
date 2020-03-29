# Merge tpms into one table

# Assign the snakemake wildcards to variables
HSrunNumberi <- snakemake@wildcards$HSrunNumber
database <- snakemake@wildcards$database
databaseRelease <- snakemake@wildcards$databaseRelease

# Load dplyr package
library(dplyr)

# Import samples information
samples <- read.delim('samples/samples.tsv')

# Select the samples for the HS run of interest
samplesRun <- dplyr::filter(samples, HSrunNumber == HSrunNumberi)

# Create a variable with sample names
sampleNames <- paste0(samplesRun$stage, '_', samplesRun$bacterium, '_', samplesRun$replicate)

# Put the folder in which the tpms are stored into a variable
folder <- paste0('output/HS', HSrunNumberi, '/counts/', database, '/release', databaseRelease, '/kallisto/')

wormStrain1 <- samplesRun[1, 'wormStrain']
stage1 <- samplesRun[1, 'stage']
bacterium1 <- samplesRun[1, 'bacterium']
replicate1 <- samplesRun[1, 'replicate']

# Store the path to the tpms of the first sample in a variable
path2abundances1 <- paste0(folder, wormStrain1, '/', stage1, '/', bacterium1, '/repl', replicate1, '/abundance.tsv')
# Read in the tpms of sample 1
abundances <- read.delim(path2abundances1)
# Rename the tpm column to the sample name
colnames(abundances)[5] <- sampleNames[1]

# Loop through all the tpm tables and add the tpm
# column to the 'abundances' data frame
for (i in 2:nrow(samplesRun)) {
  wormStraini <- samplesRun[i, 'wormStrain']
  stagei <- samplesRun[i, 'stage']
  bacteriumi <- samplesRun[i, 'bacterium']
  replicatei <- samplesRun[i, 'replicate']
  path2abundancesi <- paste0(folder, wormStraini, '/', stagei, '/', bacteriumi, '/repl', replicatei, '/abundance.tsv')
  abundancesi <- read.delim(path2abundancesi)[, c(1, 5)]
  # Rename the 'tpm' column to the sample name
  colnames(abundancesi)[2] <- sampleNames[i]
  abundances <- merge(abundances, abundancesi, by.x = 'target_id', by.y = 'target_id')
}

write.csv(abundances, paste0(folder, 'tpms.csv'), 
          row.names = FALSE, quote = FALSE)