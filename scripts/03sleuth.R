samples <- c('hb101_1', 'hb101_2', 'hb101_3', 'BIGb446_1', 'BIGb446_2', 'BIGb446_3')
condition <- c('hb101', 'hb101', 'hb101', 'BIGb446', 'BIGb446', 'BIGb446')
path <- c(paste0('output/counts/wormbase/kallisto/adults/hb101/', c('repl1', 'repl2', 'repl3')),
          paste0('output/counts/wormbase/kallisto/adults/BIGb446/', c('repl1', 'repl2', 'repl3')))
s2c <- data.frame(sample = samples, condition, path)
s2c$path <- as.character(s2c$path)

so <- sleuth::sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

write.csv(sleuth_significant, 'output/DGE/sleuth/significant.csv')

### Create a table of correspondence between the transcript name and the gene name
system("zgrep '>' output/transcriptome/wormbase/seq/c_elegans.PRJNA13758.WS267.mRNA_transcripts.fa.gz > output/transcriptome/wormbase/seq/c_elegans.PRJNA13758.WS267.mRNA_transcripts.txt")

tx2gene <- read.table('output/transcriptome/wormbase/seq/c_elegans.PRJNA13758.WS267.mRNA_transcripts.txt')
colnames(tx2gene) <- c('transcript', 'gene')
tx2gene$transcript <- stringr::str_replace(tx2gene$transcript, '>', '')
tx2gene$gene <- stringr::str_replace(tx2gene$gene, 'gene=', '')
write.csv(tx2gene, 'output/transcriptome/wormbase/seq/tx2gene.csv')

### Add the gene name to the list of significant genes
for (i in 1:nrow(sleuth_significant)) {
  genei <- tx2gene$gene[tx2gene$transcript == sleuth_significant$target_id[i]]
  sleuth_significant$gene[i] <- genei
}

ensembl <- biomaRt::useMart('ensembl', dataset = 'celegans_gene_ensembl')
biomartFilters <- biomaRt::listFilters(ensembl)
biomartAttributes <- biomaRt::listAttributes(ensembl)
wormGenes <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position', 'gene_biotype', 'description'), mart = ensembl)
wormGenes$chromosome_name <- factor(wormGenes$chromosome_name, levels = c('I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA'))
wormGenes$length <- wormGenes$end_position - (wormGenes$start_position - 1)

### Add the biological function of the gene to the significant genes
for (i in 1:nrow(sleuth_significant)) {
  descriptioni <- wormGenes$description[wormGenes$ensembl_gene_id == sleuth_significant$gene[i]]
  if (length(descriptioni) != 1) descriptioni <- 'NOTFOUND'
  sleuth_significant$description[i] <- descriptioni
}

system("sed 's/CHROMOSOME_//g' output/genome/wormbase/annotation/c_elegans.PRJNA13758.WS267.annotations.gff2 > output/genome/wormbase/annotation/c_elegans.PRJNA13758.WS267.annotationsRenamed.gff2")

write.csv(sleuth_significant, 'output/DGE/sleuth/significantMeta.csv')

snakemake --cores 7 output/counts/ensembl/featurecounts/{adults,embryos}/{hb101,BIGb446,BIGb468}/{repl1,repl2,repl3}/counts.txt

snakemake --cores 7 output/counts/ensembl/htseq/{adults,embryos}/{hb101,BIGb446,BIGb468}/{repl1,repl2,repl3}/counts.txt

