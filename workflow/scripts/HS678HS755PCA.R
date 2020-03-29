# PCA on samples from HS678 and HS755

library(DESeq2)
library(htmlwidgets)
library(plotly)
library(stringr)

# Import the counts table
countdata <- read.table('output/HS678HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsRenamed.txt', 
                        header = TRUE, row.names = 1)

# Remove first five columns (Chr, Start, End, Strand, Length)
countdata2 <- countdata[ , 6:ncol(countdata)]

# Convert to matrix
countdata3 <- as.matrix(countdata2)
head(countdata3)

condition <- str_sub(colnames(countdata3), 1, nchar(colnames(countdata3)) - 6)

# Create vectors that contain the HS run number, the worm strain,
# the worm stage and the bacterium
HSrunNumber <- sapply(condition, function(x) strsplit(x, '_')[[1]][1])
wormStrain <- sapply(condition, function(x) strsplit(x, '_')[[1]][2])
stage <- sapply(condition, function(x) strsplit(x, '_')[[1]][3])
bacterium <- sapply(condition, function(x) strsplit(x, '_')[[1]][4])

# Subset the counts and the variables to 
# only include embryos and one HB101 condition
selectionCriterion <- stage == 'embryos' & !(bacterium == 'hb101' & HSrunNumber == 'HS678')
countdataSubset <- countdata3[, selectionCriterion]
HSrunNumberSubset <- HSrunNumber[selectionCriterion]
wormStrainSubset <- wormStrain[selectionCriterion]
stageSubset <- stage[selectionCriterion]
bacteriumSubset <- bacterium[selectionCriterion]

coldataSubset <- data.frame(row.names = colnames(countdataSubset), wormStrainSubset, stageSubset, bacteriumSubset)

dds <- DESeqDataSetFromMatrix(countData = countdataSubset, colData = coldataSubset, design = ~wormStrainSubset + bacteriumSubset)
dds
dds2 <- DESeq2::DESeq(dds)

# Plot PCA
### Tomas: rlog for Diff gene expression and vst for everything else
rld <- DESeq2::rlog(dds2, blind = TRUE)
expr <- SummarizedExperiment::assay(rld)
head(expr)
p <- DESeq2::plotPCA(rld, intgroup = c('wormStrainSubset', 'bacteriumSubset')) +
  ggplot2::ggtitle('PCA on a subset of samples') + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'))
p
plotNames <- paste0('output/HS678HS755/images/scatter/PCA678_755embryos.', c('pdf', 'png', 'html'))
ggplot2::ggsave(plotNames[1], width = 10, height = 10)
ggplot2::ggsave(plotNames[2], width = 10, height = 10)
pPlotly <- plotly::ggplotly(p)
htmlwidgets::saveWidget(pPlotly, file.path(getwd(), plotNames[3]))

# Plot PCA after variance stabilising transformation
vsd <- DESeq2::varianceStabilizingTransformation(dds2)
p <- DESeq2::plotPCA(vsd, intgroup = c('wormStrainSubset', 'bacteriumSubset')) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'))
p
plotNames <- paste0('output/HS678HS755/images/scatter/PCAafterVST678_755embryos.', c('pdf', 'png', 'html'))
ggplot2::ggsave(plotNames[1], width = 10, height = 10)

# Subset the counts and the variables to 
# only include embryos and one HB101 condition
# and no mutant wdr-23 (nob27 strain)
selectionCriterion <- stage == 'embryos' & !(bacterium == 'hb101' & HSrunNumber == 'HS678') & wormStrain != 'nob27'
countdataSubset <- countdata3[, selectionCriterion]
HSrunNumberSubset <- HSrunNumber[selectionCriterion]
wormStrainSubset <- wormStrain[selectionCriterion]
stageSubset <- stage[selectionCriterion]
bacteriumSubset <- bacterium[selectionCriterion]

coldataSubset <- data.frame(row.names = colnames(countdataSubset), stageSubset, bacteriumSubset)

dds <- DESeqDataSetFromMatrix(countData = countdataSubset, colData = coldataSubset, design = ~bacteriumSubset)
dds
dds2 <- DESeq2::DESeq(dds)

# Plot PCA
### Tomas: rlog for Diff gene expression and vst for everything else
rld <- DESeq2::rlog(dds2, blind = TRUE)
expr <- SummarizedExperiment::assay(rld)
head(expr)
p <- DESeq2::plotPCA(rld, intgroup = c('bacteriumSubset')) +
  ggplot2::ggtitle('PCA on embryos without nob27') + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'))
p
plotNames <- paste0('output/HS678HS755/images/scatter/PCA678_755embryosNoNOB27.', c('pdf', 'png', 'html'))
ggplot2::ggsave(plotNames[1], width = 10, height = 10)
ggplot2::ggsave(plotNames[2], width = 10, height = 10)
pPlotly <- plotly::ggplotly(p)
htmlwidgets::saveWidget(pPlotly, file.path(getwd(), plotNames[3]))

# Plot PCA after variance stabilising transformation
vsd <- DESeq2::varianceStabilizingTransformation(dds2)
p <- DESeq2::plotPCA(vsd, intgroup = c('bacteriumSubset')) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'))
p
plotNames <- paste0('output/HS678HS755/images/scatter/PCAafterVST678_755embryosNoNOB27.', c('pdf', 'png', 'html'))
ggplot2::ggsave(plotNames[1], width = 10, height = 10)

# p$data %>% 
#   ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, col = condition, group = condition)) + 
#   ggplot2::geom_point() + 
#   ggplot2::scale_colour_manual(name="",  values = c('black', rainbow(7))) + 
#   ggplot2::xlab('First principal component (74%)') + 
#   ggplot2::ylab('Second principal component (7%)') + 
#   ggplot2::ggtitle('Principal components analysis') + 
#   ggplot2::geom_line()
# ggplot2::ggsave('output/images/scatter/cel/ref/rna/illumina/pca/pooled/noERCC/rnaExprStarDESeq2vsd.pdf')
