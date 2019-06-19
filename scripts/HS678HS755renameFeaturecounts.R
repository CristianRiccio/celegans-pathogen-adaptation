# Rename the names of the columns of the counts table to the sample name
# (variable1_variable2_variable3_replicate)

# Assign the wildcards to variables
HSrunNumberi <- '678HS755'
# HSrunNumberi <- snakemake@wildcards$HSrunNumber
databasei <- 'ensembl'
# databasei <- snakemake@wildcards$database
databaseReleasei <- '96'
# databaseReleasei <- snakemake@wildcards$databaseRelease
strandedness <- 'unstranded'
# strandedness <- snakemake@wildcards$strandedness

# Determine what the variables are
# Read in the conditions for this HS run and select the variables
conditions <- read.delim('samples/conditions.tsv')

variables <- c()

# Loop over the condition variables and find the 2 variables in 
# this project
for (i in 2:4) {
  if (apply(conditions, 2, function(x) length(unique(x)))[i] != 1) {
    variables <- c(variables, colnames(conditions)[i])
  }
}

# Print a message saying what the variables are for this HS run
print(paste0('For HS run ', HSrunNumberi, ' the variables are "', 
             paste(variables, collapse = ' ')))

# Import the counts table
path2counts <- paste0('output/HS', HSrunNumberi,
                      '/counts/featurecounts/hisat2/celegans/',
                      databasei, '/release', databaseReleasei,
                      '/minoverlap1', strandedness,
                      'FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt')
countdata <- read.table(path2counts, header = TRUE)
# countdata <- read.table(snakemake@input$counts, header = TRUE)


# Loop through the sample columns and rename the column name
for (i in 7:ncol(countdata)) {
  
  sampleCountsi <- stringr::str_split(colnames(countdata), '[.]')[[i]]
  
  wormStraini <- sampleCountsi[10]
  stagei <- sampleCountsi[11]
  bacteriumi <- sampleCountsi[12]
  
  # Assign the replicate for this sample
  replicatei <- sampleCountsi[13]
  
  # New column name (= condition name)
  conditioni <- paste0(wormStraini, '_', stagei, '_', bacteriumi, '_', replicatei)
  
  colnames(countdata)[i] <- conditioni
}

# Export the counts table with new column names
path2countsRenamed <- paste0('output/HS', HSrunNumberi,
                      '/counts/featurecounts/hisat2/celegans/',
                      databasei, '/release', databaseReleasei,
                      '/minoverlap1', strandedness,
                      'FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsRenamed.txt')
countdata <- write.table(countdata, path2countsRenamed, row.names = FALSE, quote = FALSE)
# countdata <- write.table(countdata, snakemake@output$counts, row.names = FALSE, quote = FALSE)
