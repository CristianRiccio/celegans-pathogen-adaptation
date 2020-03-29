# Keti's GO analysis

# Load R packages
library(ALL)
library(AnnotationDbi)
library(GO.db)
library(Rgraphviz)
library(topGO)

# Load the datasets
data(ALL)
data(geneList)

# Load the Gene ID/GO term map object
load('output/GO/mmusculus_gene_ensembl/idGoMap.rda')
head(idGoMap)

# Read in the gene metadata table
geneTable <- read.csv('output/GO/mmusculus_gene_ensembl/geneTable.csv')

pvalues <- read.csv(file = 'input/keti/pvalues.csv')
colnames(pvalues)[1] <- 'uniprotID'

proteinList <- merge(pvalues, geneTable, by.x = 'uniprotID', by.y = 'uniprotswissprot')

geneList <- proteinList$adj.P.Val
names(geneList) <- proteinList$ensembl_gene_id
geneList[is.na(geneList)] <- 1

sampleGOdata <- new('topGOdata', description = 'Simple session', ontology = 'BP', allGenes = geneList, geneSel = topDiffGenes, nodeSize = 1, annot = annFUN.gene2GO, gene2GO = idGoMap)

resultKS <- runTest(sampleGOdata, algorithm = 'classic', statistic = 'ks')
resultKS.elim <- runTest(sampleGOdata, algorithm = 'elim', statistic = 'ks')

pdf('output/images/graph/keti.pdf')
topGO::showSigOfNodes(sampleGOdata, score(resultKS), firstSigNodes = 5, useInfo = 'all')
# title('Graph indicating overrepresented GO terms for Keti')
dev.off()
