# Rename the names of the columns of the counts table to the sample name
# (variable1_variable2_replicate)

# Load the packages
library(dplyr)
library(stringr)

# Assign the wildcards to variables
# HSrunNumberi <- '678'
# HSrunNumberi <- '678HS755'
HSrunNumberi <- snakemake@wildcards$HSrunNumber
# databasei <- 'ensembl'
databasei <- snakemake@wildcards$database
# databaseReleasei <- '96'
databaseReleasei <- snakemake@wildcards$databaseRelease
# strandedness <- 'unstranded'
strandedness <- snakemake@wildcards$strandedness

# Read in the conditions for this HS run
conditions <- read.delim('samples/conditions.tsv')
conditionsHSrun <- conditions
if (!is.na(as.integer(HSrunNumberi))) {
  conditionsHSrun <- dplyr::filter(conditions, HSrunNumber == HSrunNumberi)
}

# Determine what the variables are
variables <- c()

# Loop over the condition variables and find the 2 variables in 
# this project
for (i in 1:4) {
  if (apply(conditionsHSrun, 2, function(x) length(unique(x)))[i] != 1) {
    variables <- c(variables, colnames(conditionsHSrun)[i])
  }
}

# Print a message saying what the variables are for this HS run
print(paste0('For HS run ', HSrunNumberi, ' the variables are "', 
             paste0(variables, collapse = ' ')))

# Import the counts table
# path2counts <- paste0('output/HS', HSrunNumberi,
#                       '/counts/featurecounts/hisat2/celegans/',
#                       databasei, '/release', databaseReleasei,
#                       '/minoverlap1', strandedness,
#                       'FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt')
# countdata <- read.table(path2counts, header = TRUE)
countdata <- read.table(snakemake@input$counts, header = TRUE)


# Loop through the sample columns and rename the column name
for (i in 7:ncol(countdata)) {
  print(paste0('column ', i))
  sampleCountsi <- stringr::str_split(colnames(countdata), '[.]')[[i]]
  variablesi <- character(length(variables))
  for (j in 1:length(variables)) {
    
    # Assign variable 1
    print(paste0('variable ', j))
    if ('HSrunNumber' == variables[j]) {
      variablesi[j] <- sampleCountsi[2]
    }
    if ('wormStrain' == variables[j]) {
      variablesi[j] <- sampleCountsi[10]
    }
    if ('stage' == variables[j]) {
      variablesi[j] <- sampleCountsi[11]
    }
    if ('bacterium' == variables[j]) {
      variablesi[j] <- sampleCountsi[12]
    }
  }
  # Assign the replicate for this sample
  replicatei <- sampleCountsi[13]
  
  # New column name (= condition name)
  conditioni <- paste(c(variablesi, replicatei), collapse = '_')
  
  colnames(countdata)[i] <- conditioni
}

# Export the counts table with new column names
# path2countsRenamed <- paste0('output/HS', HSrunNumberi,
#                       '/counts/featurecounts/hisat2/celegans/',
#                       databasei, '/release', databaseReleasei,
#                       '/minoverlap1', strandedness,
#                       'FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsRenamed.txt')
# countdata <- write.table(countdata, path2countsRenamed, row.names = FALSE, quote = FALSE)
countdata <- write.table(countdata, snakemake@output$counts, row.names = FALSE, quote = FALSE)
