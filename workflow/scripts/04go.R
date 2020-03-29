# GO term analysis

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
load('output/robjects/idGoMap.rda')
head(idGoMap)

# Adults

# Load the differentially expressed genes compared to HB101
load('output/robjects/adultsDGE.rda')

# http://www.rpubs.com/aemoore62/TopGo_colMap_Func_Troubleshoot
# solves the 'colMap' not found problem
#Defined colMap, ref. https://github.com/Bioconductor-mirror/topGO/blob/master/vignettes/topGO.Rnw
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

### Create graphs
for (pathogen in c('BIGb446adults', 'BIGb468adults')) {
  geneList <- dge[[pathogen]]$padj
  names(geneList) <- rownames(dge[[pathogen]])
  
  #  a selection function used for defining the list of differentially
  # expressed genes is also loaded under the name of
  # topDiffGenes
  # .  The function assumes that the provided
  # argument is a named vector of
  # p-values
  
  
  geneList[is.na(geneList)] <- 1
  sampleGOdata <- new('topGOdata', description = 'Simple session', ontology = 'BP', allGenes = geneList, geneSel = topDiffGenes, nodeSize = 1, annot = annFUN.gene2GO, gene2GO = idGoMap)
  
  resultKS <- runTest(sampleGOdata, algorithm = 'classic', statistic = 'ks')
  resultKS.elim <- runTest(sampleGOdata, algorithm = 'elim', statistic = 'ks')
  # pValue.classic <- score(resultKS)
  # pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
  # gstat <- termStat(sampleGOdata, names(pValue.classic))
  # gSize <- gstat$Annotated / max(gstat$Annotated) * 4
  
  
  
  # 
  #   gCol <- colMap(gstat$Significant)
  #   plot(pValue.classic, pValue.elim, xlab = 'p-value classic', ylab = 'p-value   elim', pch = 19, cex = gSize, col = gCol)
  
  # sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
  # cbind(termStat(sampleGOdata, sel.go),
  # elim = pValue.elim[sel.go],
  # classic = pValue.classic[sel.go])
  
  pdf(paste0('output/images/graph/',   pathogen, '.pdf'))
  topGO::showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5,   useInfo = 'all')
  title(paste0('Graph indicating overrepresented GO terms for ', pathogen))
  dev.off()
}

#### Embryos ####

### Load the differentially expressed genes compared to HB101
load('output/robjects/embryosDGE.rda')

# http://www.rpubs.com/aemoore62/TopGo_colMap_Func_Troubleshoot
# solves the 'colMap' not found problem
#Defined colMap, ref. https://github.com/Bioconductor-mirror/topGO/blob/master/vignettes/topGO.Rnw
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

### Create graphs
for (pathogen in c('BIGb446embryos', 'BIGb468embryos')) {
  geneList <- dge[[pathogen]]$padj
  names(geneList) <- rownames(dge[[pathogen]])
  
  #  a selection function used for defining the list of differentially
  # expressed genes is also loaded under the name of
  # topDiffGenes
  # .  The function assumes that the provided
  # argument is a named vector of
  # p-values
  
  
  geneList[is.na(geneList)] <- 1
  sampleGOdata <- new('topGOdata', description = 'Simple session', ontology = 'BP', allGenes = geneList, geneSel = topDiffGenes, nodeSize = 1, annot = annFUN.gene2GO, gene2GO = idGoMap)
  
  resultKS <- runTest(sampleGOdata, algorithm = 'classic', statistic = 'ks')
  resultKS.elim <- runTest(sampleGOdata, algorithm = 'elim', statistic = 'ks')
  # pValue.classic <- score(resultKS)
  # pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
  # gstat <- termStat(sampleGOdata, names(pValue.classic))
  # gSize <- gstat$Annotated / max(gstat$Annotated) * 4
  
  
  
  # 
  #   gCol <- colMap(gstat$Significant)
  #   plot(pValue.classic, pValue.elim, xlab = 'p-value classic', ylab = 'p-value   elim', pch = 19, cex = gSize, col = gCol)
  
  # sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
  # cbind(termStat(sampleGOdata, sel.go),
  # elim = pValue.elim[sel.go],
  # classic = pValue.classic[sel.go])
  
  pdf(paste0('output/images/graph/',   pathogen, 'embryos.pdf'))
  topGO::showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5,   useInfo = 'all')
  title(paste0('Graph indicating overrepresented GO terms for ', pathogen, ' in embryos'))
  dev.off()
}
