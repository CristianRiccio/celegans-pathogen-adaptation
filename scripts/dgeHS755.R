# DGE HS755

# Load packages
library(DESeq2)
library(purrr)
library(stringr)

source('scripts/01createFunctions.R')

# Import the counts table
countdata <- read.table('output/HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt', 
                        header = TRUE, row.names = 1)

# Remove first five columns (Chr, Start, End, Strand, Length)
countdata2 <- countdata[ , 6:ncol(countdata)]

stringr::str_split(colnames(countdata), '.')

# Remove .bam or .sam from filenames
colnames(countdata2) <- gsub('\\.[sb]am$', '', colnames(countdata2))

# Convert to matrix
countdata3 <- as.matrix(countdata2)
head(countdata3)

condition <- rep(c('n2_hb101', 'n2_pa14', 'n2_plum', 'nob27_hb101'), each = 3)

colnames(countdata3) <- paste0(condition, c('_1', '_2', '_3'))

coldata <- data.frame(row.names = colnames(countdata3), condition)

dds <- DESeqDataSetFromMatrix(countData = countdata3, colData = coldata, design = ~condition)
dds
dds2 <- DESeq2::DESeq(dds)

res <- DESeq2::results(dds2)
print(head(res))
res$padj %>% table %>% head

SummarizedExperiment::mcols(res, use.names = TRUE)

dge <- DESeq2::results(dds2, c('condition', 'hb101', 'plum')) %>% data.frame()

# Add the protein description as a column to the data frames

dge$description <- purrr::map_chr(rownames(dge), idToDescription)
# Add the gene name as a column to the data frame
dge$external_gene_name <- purrr::map_chr(rownames(dge), idToGeneName)
# Select significant genes (p-value < 5%)
significantGenes <- dge[which(dge$padj < 0.05), ]
print(nrow(significantGenes))
# Write out the list of significant genes
write.csv(significantGenes, paste0('output/HS755/DGE/hb101plum.csv'))
# Save the DGE analysis in an R object
save(dge, file = 'output/HS755/robjects/hb101plum.rda')

# Plot PCA
### Tomas: rlog for Diff gene expression and vst for everything else
rld <- DESeq2::rlog(dds2, blind = TRUE)
expr <- SummarizedExperiment::assay(rld)
head(expr)
DESeq2::plotPCA(rld, intgroup = ('condition'))
ggplot2::ggsave('output/HS755/images/scatter/PCA755.pdf',
                width = 10, height = 10)

