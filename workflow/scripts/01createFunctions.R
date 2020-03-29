# Create functions that convert gene IDs into gene names and descriptions

# Import packages
library(biomaRt)
library(tidyr)

# biomartDataset <- 'celegans_gene_ensembl'
# biomartDataset <- 'mmusculus_gene_ensembl'
biomartDataset <- snakemake@wildcards$biomartDataset

ensembl <- biomaRt::useMart('ensembl', dataset = biomartDataset)
biomartFilters <- biomaRt::listFilters(ensembl)
biomartAttributes <- biomaRt::listAttributes(ensembl)
wormGenes <- biomaRt::getBM(attributes = c('external_gene_name', 'chromosome_name', 'gene_biotype', 'description', 'ensembl_gene_id', 'go_id', 'name_1006', 'start_position', 'end_position', 'uniprotswissprot'), mart = ensembl)
wormGenes$chromosome_name <- factor(wormGenes$chromosome_name, levels = c('I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA'))
save(wormGenes, file = paste0('output/GO/', biomartDataset, '/table.rda'))

wormGenes$length <- wormGenes$end_position - (wormGenes$start_position - 1)

idToDescription <- function(ensemblGeneId = 'WBGene00000001') {
  return(wormGenes$description[wormGenes$ensembl_gene_id == ensemblGeneId][1])
}

idToGeneName <- function(ensemblGeneId = 'WBGene00000001') {
  return(wormGenes$external_gene_name[wormGenes$ensembl_gene_id == ensemblGeneId][1])
}