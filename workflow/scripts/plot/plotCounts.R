# Plot featurecounts

# Load R packages
library(ggplot2)
library(htmlwidgets)
library(plotly)
library(stringr)
library(tidyr)

# Making the plots interactive is a nice idea but
# they require a lot of computational resources to 
# interact with, so limit interaction to plots 
# with little data

# Assign the wildcards to variables
# HSrunNumberi <- 678
HSrunNumberi <- snakemake@wildcards$HSrunNumber

# Import the counts table
# path2counts <- paste0('output/HS', HSrunNumberi, '/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsRenamed.txt')
path2counts <- snakemake@input$counts
counts <- read.table(path2counts, header = TRUE)

# Tidy the counts so that each row represent a count in a certain sample
countsTidy <- tidyr::gather(counts, 'sample', 'counts', -Geneid, -Chr, -Start, -End, -Strand, -Length)

# Store the sample variables (2 variables + replicates) in an R variable 
sampleVar <- stringr::str_split(countsTidy$sample, '_', simplify = TRUE)

# Add columns to the data frame: 2 variables + replicate
# Determine what the variables are
# Read in the conditions for this HS run and select the two variables
# conditions <- read.delim('samples/conditions.tsv')
conditions <- read.delim(snakemake@input$conditions)
conditionsHSrun <- dplyr::filter(conditions, HSrunNumber == HSrunNumberi)

variables <- c()

# Loop over the condition variables and find the 2 variables in 
# this project
for (i in 2:4) {
  if (apply(conditionsHSrun, 2, function(x) length(unique(x)))[i] != 1) {
    variables <- c(variables, colnames(conditionsHSrun)[i])
  }
}

countsTidy[[variables[1]]] <- sampleVar[, 1]
countsTidy[[variables[2]]] <- sampleVar[, 2]
countsTidy$replicate <- sampleVar[, 3]

# Summarise the counts and export them to a file
# path2summary <- paste0('output/HS', HSrunNumberi, '/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsSummary.txt')
path2summary <- snakemake@output$summary
sink(path2summary)
countsSummary <- summary(countsTidy$counts)
names(countsSummary) <- c('min', 'quant1', 'median', 'mean', 'quant3', 'max')
countsSummary
sink()

# Calculate the total number of counts in each sample
countsPerSample <- dplyr::group_by(countsTidy, .dots = c(variables[1], variables[2], 'replicate')) %>% 
  dplyr::summarise(totalCounts = sum(counts)) 

# Write out a table of the total number of counts in each sample
# path2table <- paste0('output/HS', HSrunNumberi, '/tables/celegans/totalCounts.csv')
path2table <- snakemake@output$countsTable
write.csv(countsPerSample, path2table, row.names = FALSE)

# Plot the total number of counts in each sample
p <- ggplot2::ggplot(countsPerSample, ggplot2::aes(x = replicate, y = totalCounts)) +
  ggplot2::geom_bar(stat = 'identity') + 
  ggplot2::ylab('Total counts / millions') +
  ggplot2::ggtitle('Total counts in each sample') +
  ggplot2::facet_grid(reformulate(variables[1], variables[2])) +
  ggplot2::scale_y_continuous(labels = function(x) x/1e6)
p
# plotNames <- paste0('output/HS', HSrunNumberi, '/images/bar/celegans/totalCounts.', c('pdf', 'png', 'html'))
print(snakemake@output$totalCounts)
plotNames <- snakemake@output$totalCounts
ggplot2::ggsave(plotNames[1], height = 6, width = 9)
ggplot2::ggsave(plotNames[2], height = 6, width = 9)
pPlotly <- plotly::ggplotly(p)
pPlotly
htmlwidgets::saveWidget(pPlotly, file.path(getwd(), plotNames[3]))


# Plot the distribution of counts
# Boxplots
p <- ggplot2::ggplot(countsTidy, ggplot2::aes(x = replicate, y = counts)) +
  ggplot2::geom_boxplot() +
  ggplot2::ylab('Counts') +
  ggplot2::ggtitle('Number of counts across all genes and samples') +
  ggplot2::facet_grid(reformulate(variables[1], variables[2]))
p
# plotNames <- paste0('output/HS', HSrunNumberi, '/images/boxplot/counts.', c('pdf', 'png', 'html'))
plotNames <- snakemake@output$boxplotCounts
ggplot2::ggsave(plotNames[1], height = 8, width = 15)
ggplot2::ggsave(plotNames[2], height = 8, width = 15)
pPlotly <- plotly::ggplotly(p)
htmlwidgets::saveWidget(pPlotly, file.path(getwd(), plotNames[3]))
# Plot the counts on a logarithmic scale
p <- p + ggplot2::scale_y_log10()
p
# plotNames <- paste0('output/HS', HSrunNumberi, '/images/boxplot/countsLog10.', c('pdf', 'png', 'html'))
plotNames <- snakemake@output$boxplotCountsLog10
ggplot2::ggsave(plotNames[1], height = 8, width = 15)
ggplot2::ggsave(plotNames[2], height = 8, width = 15)
pPlotly <- plotly::ggplotly(p)
htmlwidgets::saveWidget(pPlotly, file.path(getwd(), plotNames[3]))

# Histograms
p <- ggplot2::ggplot(countsTidy, ggplot2::aes(x = counts)) +
  ggplot2::geom_histogram() +
  ggplot2::xlab('Featurecounts') +
  ggplot2::ylab('Counts in each bin') +
  ggplot2::ggtitle('Number of counts across all genes and samples') +
  ggplot2::facet_grid(reformulate(variables[1], variables[2]))
p
# plotNames <- paste0('output/HS', HSrunNumberi, '/images/histogram/counts.', c('pdf', 'png', 'html'))
plotNames <- snakemake@output$histogramCounts
ggplot2::ggsave(plotNames[1], height = 8, width = 15)
ggplot2::ggsave(plotNames[2], height = 8, width = 15)
# Plot the counts on a logarithmic scale
p <- p + ggplot2::scale_x_log10()
p
pPlotly <- plotly::ggplotly(p)
htmlwidgets::saveWidget(pPlotly, file.path(getwd(), plotNames[3]))
# plotNames <- paste0('output/HS', HSrunNumberi, '/images/histogram/countsLog10.', c('pdf', 'png', 'html'))
plotNames <- snakemake@output$histogramCountsLog10
ggplot2::ggsave(plotNames[1], height = 8, width = 15)
ggplot2::ggsave(plotNames[2], height = 8, width = 15)
pPlotly <- plotly::ggplotly(p)
htmlwidgets::saveWidget(pPlotly, file.path(getwd(), plotNames[3]))