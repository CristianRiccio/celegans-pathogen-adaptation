# Create functions that convert gene IDs into gene names and descriptions

# Import packages
library(biomaRt)
library(tidyr)

# Create the folder where R objects will be stored
dir.create('output/robjects', recursive = TRUE)

ensembl <- biomaRt::useMart('ensembl', dataset = 'celegans_gene_ensembl')
biomartFilters <- biomaRt::listFilters(ensembl)
biomartAttributes <- biomaRt::listAttributes(ensembl)
wormGenes <- biomaRt::getBM(attributes = c('external_gene_name', 'chromosome_name', 'gene_biotype', 'description', 'ensembl_gene_id', 'go_id', 'name_1006', 'uniprot_gn', 'start_position', 'end_position'), mart = ensembl)
wormGenes$chromosome_name <- factor(wormGenes$chromosome_name, levels = c('I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA'))
save(wormGenes, file = 'output/robjects/wormGenes.rda')

wormGenes$length <- wormGenes$end_position - (wormGenes$start_position - 1)

head(wormGenes)

biomartAttributes

idToDescription <- function(ensemblGeneId = 'WBGene00000001') {
  return(wormGenes$description[wormGenes$ensembl_gene_id == ensemblGeneId][1])
}

idToGeneName <- function(ensemblGeneId = 'WBGene00000001') {
  return(wormGenes$external_gene_name[wormGenes$ensembl_gene_id == ensemblGeneId][1])
}

idToDescription()

purrr::map_chr(c('WBGene00000001', 'WBGene00000002', 'WBGene00000003'), idToDescription)


### Create an Ensembl ID <-> GO term map

idGoTable <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'go_id'), mart = ensembl)
idGoMap <- list()
for (ensemblGeneId in unique(idGoTable$ensembl_gene_id)) {
  idGoMap[[ensemblGeneId]] <- idGoTable$go_id[idGoTable$ensembl_gene_id == ensemblGeneId]
}
save(idGoTable, file = 'output/robjects/idGoTable.rda')
save(idGoMap, file = 'output/robjects/idGoMap.rda')


