### DGE analysis
setwd('~/Documents/phd/projects/nick')
library(DESeq2)

# paths <- c(paste0('output/counts/ensembl/htseq/adults/hb101/', c('repl1', 'repl2', 'repl3'), '/counts.txt'),
#            paste0('output/counts/ensembl/htseq/adults/BIGb446/', c('repl1', 'repl2', 'repl3'), '/counts.txt'))
# 
# # DESeq2::DESeqDataSetFromMatrix(countData=counts, colData=sampleInfo, design=~condition)
# 
# sampleData <- data.frame(SampleName = c('adult1', 'adult2', 'adult3', 'adult4', 'adult5', 'adult6'), 
#                          FileName = paths,
#                          Condition = c('hb101', 'hb101', 'hb101', 'BIGb446', 'BIGb446', 'BIGb446'))
# 
# dds <- DESeq2::DESeqDataSetFromHTSeqCount(sampleTable = sampleData, 
#                                           design = ~Condition)
# print(dds)
# dds$Condition <- relevel(dds$Condition, 'N2')
# # DESeq2::plotCounts(dds, 'WBGene00000001', intgroup = 'roxel')
# dds2 <- DESeq2::DESeq(dds)
# 


countdata <- read.table('output/counts/ensembl/featurecounts/counts.txt', 
                        header = TRUE, row.names = 1)

# Remove first five columns (chr, start, end, strand, length)
countdata2 <- countdata[ ,6:ncol(countdata)]

# Remove .bam or .sam from filenames
colnames(countdata2) <- gsub("\\.[sb]am$", "", colnames(countdata2))

# Convert to matrix
countdata3 <- as.matrix(countdata2)
head(countdata3)

### Let's start a DGE analysis with just 3 conditions, the adults
countdata4 <- countdata3[, 1:9]

conditions1 <- stringr::str_extract(colnames(countdata4), '[.][a-zA-Z0-9]*[.]repl')
conditions2 <- substr(conditions1, 2, nchar(conditions1) - 5)
conditions3 <- stringr::str_extract(colnames(countdata4), 'ensembl[.][a-z]*[.]')
conditions4 <- substr(conditions3, 9, nchar(conditions3) - 1)
condition <- paste0(conditions2, conditions4)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
colnames(countdata4) <- paste0(condition, c('_1', '_2', '_3'))

coldata <- data.frame(row.names=colnames(countdata4), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata4, colData=coldata, design=~condition)
dds
dds2 <- DESeq2::DESeq(dds)

### Tomas: rlog for Diff gene expression and vst for everything else
rld <- DESeq2::rlog(dds2, blind = TRUE)
expr <- SummarizedExperiment::assay(rld)
head(expr)

# install.packages('pheatmat')
library(pheatmap)

pdf('output/images/heatmap/adultCorrelations.pdf')
pheatmap::pheatmap(cor(expr), clustering_method = 'average',
                   main = 'Correlation between normalised gene expression values between samples')
dev.off()


pheatmap::pheatmap(expr, clustering_method = 'average',
                   main = 'Correlation between normalised gene expression values between samples')


### Hierarchical clustering
clusters <- hclust(dist(t(expr)))


### Plot the tree of the hierarchical clustering
pdf('output/images/tree/adultsHierarchicalClustering.pdf')
plot(clusters, main = 'Hierarchical clustering of samples 
     based on normalised gene expression values')
dev.off()

# snakemake --dag output/counts/ensembl/htseq/{adults,embryos}/{hb101,BIGb446,BIGb468}/{repl1,repl2,repl3}/counts.txt | dot -Tsvg > dag.svg
# convert dag.svg dag.pdf

library(dplyr)
source('code/01createFunctions.R')

DESeq2::plotPCA(rld, intgroup = ('condition'))
ggplot2::ggsave('output/images/scatter/adultsPCArnaExprStarDESeq2.pdf',
                width = 10, height = 10)

pca <- prcomp(t(expr))
loadings <- as.data.frame(pca$rotation)
loading1 <- loadings['PC1']
loading1ord <- loading1[order(abs(loading1$PC1), decreasing = TRUE), , drop = FALSE]
sapply(rownames(loading1ord)[1:50], idToDescription) -> pc1genes
write.csv(pc1genes, file = 'output/DGE/adultsPC1.csv')

loading2 <- loadings['PC2']
loading2ord <- loading2[order(abs(loading2$PC2), decreasing = TRUE), , drop = FALSE]
sapply(rownames(loading2ord)[1:50], idToDescription) -> pc2genes
write.csv(pc2genes, file = 'output/DGE/adultsPC2.csv')

vsd <- DESeq2::varianceStabilizingTransformation(dds2)
DESeq2::plotPCA(vsd, intgroup = ('condition')) +
  ggplot2::ggtitle('PCA of gene expression in adults after variance stabilising transformation')
ggplot2::ggsave('output/images/scatter/adultsRnaExprStarDESeq2vsd.pdf',
                width = 20, height = 20)


res <- DESeq2::results(dds2)
print(head(res))
res$padj %>% table %>% head

SummarizedExperiment::mcols(res, use.names = TRUE)

dge <- list()
significantGenes <- list()
for (group in unique(condition[condition != 'hb101adults'])) {
  dge[[group]] <- DESeq2::results(dds2, c('condition', group, 'hb101adults')) %>% as.data.frame()
  rownames(dge[[group]]) <- rownames(DESeq2::results(dds2, c('condition', group, 'hb101adults')))
  
  ### Add the protein description as a column to the data frames
  dge[[group]]$description <- purrr::map_chr(rownames(dge[[group]]), idToDescription)
  ### Add the gene name as a column to the data frame
  dge[[group]]$external_gene_name <- purrr::map_chr(rownames(dge[[group]]), idToGeneName)
  
  significantGenes[[group]] <- dge[[group]][which(dge[[group]]['padj'] < 0.05),]
  print(group)
  print(nrow(significantGenes[[group]]))
  write.csv(significantGenes[[group]], paste0('output/DGE/', group, '.csv'))
}
save(dge, file = 'output/robjects/adultsDGE.rda')

### Plot the number of significantly differentially expressed genes 
sapply(significantGenes, nrow) %>% as.data.frame() -> sigGenes
names(sigGenes) <- 'numberDEgenes'
sigGenes$strain <- rownames(sigGenes)

### Plot the number of differentially expressed genes
sigGenes  %>% 
  ggplot2::ggplot(ggplot2::aes(x = strain, y = numberDEgenes)) + 
  ggplot2::geom_bar(stat = 'identity') + ggplot2::xlab('Strain') + 
  ggplot2::ylab('Number of differentially expressed genes') +
  ggplot2::ggtitle('Number of genes differentially expressed compared to HB101') +
  ggplot2::geom_text(ggplot2::aes(label = numberDEgenes))
ggplot2::ggsave('output/images/bar/adultsNbDiffExprGenes.pdf', width = 9)


### Compile a list of gene IDs that are differentially expressed in a strain compared to HB101
mart <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
                         dataset = 'celegans_gene_ensembl',
                         host = 'ensembl.org')

biomaRt::listAttributes(mart)


t2g <- biomaRt::getBM(attributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name', 'wikigene_description'), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)


t2g$wikigene_description[t2g$ensembl_gene_id %in% c('WBGene00000001', 'WBGene00000002')]
deGenes <- c()
for (group in unique(condition[condition != 'hb101adults'])) {
  deGenes <- c(deGenes, rownames(significantGenes[[group]])
  )
}
deGenes <- unique(deGenes)

