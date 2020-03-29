# Create Ensembl ID <-> GO term map

# Import packages
library(biomaRt)

# biomartDataset <- 'celegans_gene_ensembl'
# biomartDataset <- 'mmusculus_gene_ensembl'
biomartDataset <- snakemake@wildcards$biomartDataset

ensembl <- biomaRt::useMart('ensembl', dataset = biomartDataset)

geneTable <- biomaRt::getBM(attributes = c('external_gene_name', 'chromosome_name', 'gene_biotype', 'description', 'ensembl_gene_id', 'go_id', 'name_1006', 'start_position', 'end_position', 'uniprotswissprot'), mart = ensembl)
write.csv(geneTable, file = paste0('output/GO/', biomartDataset, '/geneTable.csv'))

idGoTable <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'go_id'), mart = ensembl)
idGoMap <- list()
for (ensemblGeneId in unique(idGoTable$ensembl_gene_id)) {
  idGoMap[[ensemblGeneId]] <- idGoTable$go_id[idGoTable$ensembl_gene_id == ensemblGeneId]
}

write.csv(idGoTable, file = paste0('output/GO/', biomartDataset, '/ensembl_gene_id_go_id_table.csv'))
save(idGoMap, file = paste0('output/GO/', biomartDataset, '/idGoMap.rda'))
