### DGE embryos

# Load packages
library(DESeq2)
library(stringr)

#source('scripts/01createFunctions.R')
# paths <- c(paste0('output/counts/ensembl/htseq/embryos/hb101/', c('repl1', 'repl2', 'repl3'), '/counts.txt'),
#            paste0('output/counts/ensembl/htseq/embryos/BIGb446/', c('repl1', 'repl2', 'repl3'), '/counts.txt'))
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


countdata <- read.table(snakemake@input$counts, 
                        header = TRUE, row.names = 1)
# countdata <- read.table('output/counts/featurecounts/ensembl/release95/counts.txt', 
#                        header = TRUE, row.names = 1)
# Remove first five columns (chr, start, end, strand, length)
countdata2 <- countdata[ ,6:ncol(countdata)]

# Remove .bam or .sam from filenames
colnames(countdata2) <- gsub("\\.[sb]am$", "", colnames(countdata2))

# Convert to matrix
countdata3 <- as.matrix(countdata2)
head(countdata3)

### Let's start a DGE analysis with just 3 conditions, the embryos
countdata4 <- countdata3[, 10:18]

conditions1 <- stringr::str_extract(colnames(countdata4), '[.][a-zA-Z0-9]*[.]repl')
conditions2 <- substr(conditions1, 2, nchar(conditions1) - 5)
conditions3 <- stringr::str_extract(colnames(countdata4), 'ensembl[.][a-z]*[.]')
conditions4 <- substr(conditions3, 9, nchar(conditions3) - 1)
condition <- paste0(conditions2, conditions4)
print(3)
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
colnames(countdata4) <- paste0(condition, c('_1', '_2', '_3'))

coldata <- data.frame(row.names = colnames(countdata4), condition)
dds <- DESeqDataSetFromMatrix(countData = countdata4, colData = coldata, design = ~condition)
dds
dds2 <- DESeq2::DESeq(dds)
# save(dds2, file = 'output/robjects/featurecounts/ensembl/release95/dds2.Rda')
save(dds2, snakemake@output$dds2)

#res <- DESeq2::results(dds2)
#print(head(res))
#res$padj %>% table %>% head
#
#SummarizedExperiment::mcols(res, use.names = TRUE)
#
#dge <- list()
#significantGenes <- list()
#for (group in unique(condition[condition != 'hb101embryos'])) {
#  dge[[group]] <- DESeq2::results(dds2, c('condition', group, 'hb101embryos')) %>% as.data.frame()
#  rownames(dge[[group]]) <- rownames(DESeq2::results(dds2, c('condition', group, 'hb101embryos')))
#  
#  ### Add the protein description as a column to the data frames
#  dge[[group]]$description <- purrr::map_chr(rownames(dge[[group]]), idToDescription)
#  ### Add the gene name as a column to the data frame
#  dge[[group]]$external_gene_name <- purrr::map_chr(rownames(dge[[group]]), idToGeneName)
#  
#  significantGenes[[group]] <- dge[[group]][which(dge[[group]]['padj'] < 0.05),]
#  print(group)
#  print(nrow(significantGenes[[group]]))
#  write.csv(significantGenes[[group]], paste0('output/DGE/', group, '.csv'))
#}
#save(dge, file = 'output/robjects/embryosDGE.rda')
#
#### Plot the number of significantly differentially expressed genes 
#sapply(significantGenes, nrow) %>% as.data.frame() -> sigGenes
#names(sigGenes) <- 'numberDEgenes'
#sigGenes$strain <- rownames(sigGenes)
#
#### Plot the number of differentially expressed genes
#sigGenes  %>% 
#  ggplot2::ggplot(ggplot2::aes(x = strain, y = numberDEgenes)) + 
#  ggplot2::geom_bar(stat = 'identity') + ggplot2::xlab('Strain') + 
#  ggplot2::ylab('Number of differentially expressed genes') +
#  ggplot2::ggtitle('Number of genes differentially expressed compared to HB101') +
#  ggplot2::geom_text(ggplot2::aes(label = numberDEgenes))
#ggplot2::ggsave('output/images/bar/embryosNbDiffExprGenes.pdf', width = 9)


### Compile a list of gene IDs that are differentially expressed in a strain compared to HB101
#mart <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
#                         dataset = 'celegans_gene_ensembl',
#                         host = 'ensembl.org')
#
#biomaRt::listAttributes(mart)
#
#
#t2g <- biomaRt::getBM(attributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name', 'wikigene_description'), mart = mart)
#t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
#                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
#
#
#t2g$wikigene_description[t2g$ensembl_gene_id %in% c('WBGene00000001', 'WBGene00000002')]
#deGenes <- c()
#for (group in unique(condition[condition != 'hb101embryos'])) {
#  deGenes <- c(deGenes, rownames(significantGenes[[group]])
#  )
#}
#deGenes <- unique(deGenes)

### Plot the expression values of one condition on the x-axis and the expression 
### values of another condition on the y-axis
#exprDf <- as.data.frame(expr)
#exprDf$gene <- rownames(exprDf)
#exprTidy <- tidyr::gather(exprDf, sample, expression, -gene) %>% 
#  dplyr::arrange(gene, sample)
#exprTidy$condition <- substr(exprTidy$sample, 1, nchar(exprTidy$sample) - 2)
#
#meanExpr <- dplyr::group_by(exprTidy, condition, gene) %>% 
#  dplyr::summarise(meanExpr = mean(expression)) %>% 
#  dplyr::arrange(gene)
#
#meanExprs <- tidyr::spread(meanExpr, condition, meanExpr)
#
#ggplot2::ggplot(meanExprs, ggplot2::aes(x = hb101embryos, y = BIGb446embryos)) +
#  ggplot2::geom_point()
#ggplot2::ggsave('output/images/scatter/BIGb446embryosVShb101embryos.pdf')
#
## log2foldchange vs base mean expression
#### Plot the log2FoldChange vs the basemean of genes
#### but first identify and remove the most outrageously expressed genes
#comparison <- 'BIGb446embryos_HB101embryos'
#dge$BIGb446embryos$signif <- 'NS'
#dge$BIGb446embryos$signif[dge$BIGb446embryos$padj < 0.05 & abs(dge$BIGb446embryos$log2FoldChange) > 1] <- 'signif'
#topRemoved <- 10
#dplyr::arrange(dge$BIGb446embryos, desc(baseMean))[topRemoved:nrow(dge$BIGb446embryos), ] %>%
#  ggplot2::ggplot(ggplot2::aes(x = baseMean, y = log2FoldChange, col = signif, alpha = 0.001)) +
#  ggplot2::geom_point() +
#  ggplot2::scale_color_manual(values=c("black", "red")) +
#  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
#                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
#  ggplot2::xlab('Base mean expression level') +
#  ggplot2::ylab('Log2-fold change in mean expression') +
#  ggplot2::ggtitle(paste0('Differentially gene expression in ', comparison))
#plotName <- paste0('output/images/scatter/release94/log2FoldChangeVSbaseMean', comparison, 'removedTop', topRemoved)
#ggplot2::ggsave(paste0(plotName, '.pdf'), width = 10)
#ggplot2::ggsave(paste0(plotName, '.png'), width = 10)
#ggplot2::ggsave(paste0(plotName, '.jpg'), width = 10)

# log2foldchange vs base mean expression
### Plot the log2FoldChange vs the basemean of genes
### but first identify and remove the most outrageously expressed genes
#comparison <- 'BIGb468embryos_HB101embryos'
#dge$BIGb468embryos$signif <- 'NS'
#dge$BIGb468embryos$signif[dge$BIGb468embryos$padj < 0.05 & abs(dge$BIGb468embryos$log2FoldChange) > 1] <- 'signif'
#topRemoved <- 10
#dplyr::arrange(dge$BIGb468embryos, desc(baseMean))[topRemoved:nrow(dge$BIGb468embryos), ] %>%
#  ggplot2::ggplot(ggplot2::aes(x = baseMean, y = log2FoldChange, col = signif, alpha = 0.001)) +
#  ggplot2::geom_point() +
#  ggplot2::scale_color_manual(values=c("black", "red")) +
#  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
#                 panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
#  ggplot2::xlab('Base mean expression level') +
#  ggplot2::ylab('Log2-fold change in mean expression') +
#  ggplot2::ggtitle(paste0('Differentially gene expression in ', comparison))
#plotName <- paste0('output/images/scatter/release94/log2FoldChangeVSbaseMean', comparison, 'removedTop', topRemoved)
#ggplot2::ggsave(paste0(plotName, '.pdf'), width = 10)
#ggplot2::ggsave(paste0(plotName, '.png'), width = 10)
#ggplot2::ggsave(paste0(plotName, '.jpg'), width = 10)
