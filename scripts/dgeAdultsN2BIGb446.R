# DGE N2 on HB101 and BIGb446 adults

# Load packages
library(DESeq2)
library(purrr)
library(stringr)

source('scripts/01createFunctions.R')

# Import the counts table
countdata <- read.table('output/HS678/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt', 
                        header = TRUE, row.names = 1)

# Remove first five columns (Chr, Start, End, Strand, Length)
countdata2 <- countdata[ , 6:ncol(countdata)]

# Convert to matrix
countdata3 <- as.matrix(countdata2)
head(countdata3)

# Subset the counts to include only HB101 and BIGb446
countdata4 <- countdata3[, 1:6]

condition <- rep(c('adults_hb101', 'adults_BIGb446'), each = 3)

colnames(countdata4) <- paste0(condition, c('_1', '_2', '_3'))

coldata <- data.frame(row.names = colnames(countdata4), condition)

dds <- DESeqDataSetFromMatrix(countData = countdata4, colData = coldata, design = ~condition)
dds
dds2 <- DESeq2::DESeq(dds)

res <- DESeq2::results(dds2)
print(head(res))
res$padj %>% table %>% head

SummarizedExperiment::mcols(res, use.names = TRUE)

dge <- DESeq2::results(dds2, c('condition', 'adults_hb101', 'adults_BIGb446')) %>% data.frame()

# Add the protein description as a column to the data frames

dge$description <- purrr::map_chr(rownames(dge), idToDescription)
# Add the gene name as a column to the data frame
dge$external_gene_name <- purrr::map_chr(rownames(dge), idToGeneName)
# Select significant genes (p-value < 5%)
significantGenes <- dge[which(dge$padj < 0.05), ]
print(nrow(significantGenes))
# Write out the list of significant genes
write.csv(significantGenes, paste0('output/HS678/DGE/adultsHB101BIGb446.csv'))
# Save the DGE analysis in an R object
save(dge, file = 'output/HS678/robjects/adultsHB101BIGb446.rda')

