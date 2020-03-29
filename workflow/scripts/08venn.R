### Genes differentially expressed in pathogen-exposed worms compared to HB101-exposed worms

### Load the package to draw Venn diagrams
library(VennDiagram)

### Adults
### Load the tables containing the differentially expressed genes
BIGb446adults <- read.csv('output/DGE/BIGb446adults.csv')
colnames(BIGb446adults)[1] <- 'geneID'
BIGb468adults <- read.csv('output/DGE/BIGb468adults.csv')
colnames(BIGb468adults)[1] <- 'geneID'
### Count the number of differentially expressed genes in each condition
BIGb446adultsN <- nrow(BIGb446adults)
BIGb468adultsN <- nrow(BIGb468adults)
### Count the number of genes that are differentially expressed in both conditions
geneOverlap <- intersect(BIGb446adults$geneID, BIGb468adults$geneID)


grid.newpage()

### Plot the Venn diagrams with the labels
pdf('output/images/venn/BIGb446BIGb468adultsDGElabelled.pdf')
g <- draw.pairwise.venn(BIGb446adultsN, BIGb468adultsN, length(geneOverlap),
                        category = (c('BIGb446adultsN','BIGb468adultsN')),
                        lty = rep('blank',2), fill = c('blue', 'red'),
                        alpha = c(0.5, 0.5), cex = 1, cat.fontface = 2)
grid.text('Genes differentially expressed in pathogen-exposed adults 
          compared to HB101-exposed adults', 
          y = 0.97, gp = gpar(col='black', cex=1))
dev.off()

### Plot the Venn diagrams without the labels so that they can be added
### downstream more freely
pdf('output/images/venn/BIGb446BIGb468adultsDGEunlabelled.pdf')
g <- draw.pairwise.venn(BIGb446adultsN, BIGb468adultsN, 
                        length(geneOverlap), category = (c('', '')), 
                        lty = rep('blank',2), fill = c('blue', 'red'),  
                        alpha = c(0.5, 0.5), cex = 1,cat.fontface = 2)
grid.text('', y=0.97, gp = gpar(col = 'black', cex=2))
dev.off()

### Embryos
### Load the tables containing the differentially expressed genes
BIGb446embryos <- read.csv('output/DGE/BIGb446embryos.csv')
colnames(BIGb446embryos)[1] <- 'geneID'
BIGb468embryos <- read.csv('output/DGE/BIGb468embryos.csv')
colnames(BIGb468embryos)[1] <- 'geneID'
### Count the number of differentially expressed genes in each condition
BIGb446embryosN <- nrow(BIGb446embryos)
BIGb468embryosN <- nrow(BIGb468embryos)
### Count the number of genes that are differentially expressed in both conditions
geneOverlap <- intersect(BIGb446embryos$geneID, BIGb468embryos$geneID)


grid.newpage()

### Plot the Venn diagrams with the labels
pdf('output/images/venn/BIGb446BIGb468embryosDGElabelled.pdf')
g <- draw.pairwise.venn(BIGb446embryosN, BIGb468embryosN, length(geneOverlap),
                        category = (c('BIGb446embryosN','BIGb468embryosN')),
                        lty = rep('blank',2), fill = c('blue', 'red'),
                        alpha = c(0.5, 0.5), cex = 1, cat.fontface = 2)
grid.text('Genes differentially expressed in pathogen-exposed embryos 
          compared to HB101-exposed embryos', 
          y = 0.97, gp = gpar(col='black', cex=1))
dev.off()

### Plot the Venn diagrams without the labels so that they can be added
### downstream more freely
pdf('output/images/venn/BIGb446BIGb468embryosDGEunlabelled.pdf')
g <- draw.pairwise.venn(BIGb446embryosN, BIGb468embryosN, 
                        length(geneOverlap), category = (c('', '')), 
                        lty = rep('blank',2), fill = c('blue', 'red'),  
                        alpha = c(0.5, 0.5), cex = 1,cat.fontface = 2)
grid.text('', y=0.97, gp = gpar(col = 'black', cex=2))
dev.off()
