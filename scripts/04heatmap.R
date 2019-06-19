library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(pheatmap)

### Tomas: rlog for Diff gene expression and vst for everything else
rld <- DESeq2::rlog(dds2, blind = TRUE)
expr <- SummarizedExperiment::assay(rld)

##pdf(snakemake@output$heatmap)
#pdf('output/images/heatmap/ensembl/release95/adultCorrelations.pdf')
#pheatmap::pheatmap(cor(expr), clustering_method = 'average',
#                   main = 'Correlation between normalised gene expression values between samples')
#dev.off()
#

# pheatmap::pheatmap(expr, clustering_method = 'average',
#                    main = 'Correlation between normalised gene expression values between samples')


### Hierarchical clustering
#clusters <- hclust(dist(t(expr)))
#
#
#### Plot the tree of the hierarchical clustering
#pdf('output/images/tree/embryosHierarchicalClustering.pdf')
#plot(clusters, main = 'Hierarchical clustering of samples 
#     based on normalised gene expression values')
#dev.off()
#
#
#DESeq2::plotPCA(rld, intgroup = ('condition'))
#ggplot2::ggsave('output/images/scatter/embryosPCArnaExprStarDESeq2.pdf',
#                width = 10, height = 10)
#
#pca <- prcomp(t(expr))
#loadings <- as.data.frame(pca$rotation)
#loading1 <- loadings['PC1']
#loading1ord <- loading1[order(abs(loading1$PC1), decreasing = TRUE), , drop = FALSE]
#sapply(rownames(loading1ord)[1:50], idToDescription) -> pc1genes
#write.csv(pc1genes, file = 'output/DGE/embryosPC1.csv')
#
#loading2 <- loadings['PC2']
#loading2ord <- loading2[order(abs(loading2$PC2), decreasing = TRUE), , drop = FALSE]
#sapply(rownames(loading2ord)[1:50], idToDescription) -> pc2genes
#write.csv(pc2genes, file = 'output/DGE/embryosPC2.csv')
#
#vsd <- DESeq2::varianceStabilizingTransformation(dds2)
#DESeq2::plotPCA(vsd, intgroup = ('condition')) +
#  ggplot2::ggtitle('PCA of gene expression in embryos after variance stabilising transformation')
#ggplot2::ggsave('output/images/scatter/embryosRnaExprStarDESeq2vsd.pdf',
#                width = 20, height = 20)
