tpms <- read.table('output/counts/ensembl/featurecounts/counts.txt',
                   header = TRUE, row.names = 'Geneid')

# Remove first five columns (chr, start, end, strand, length)
justTpms <- tpms[ ,6:ncol(tpms)]

### Let's start a DGE analysis with just 3 conditions, the adults
tpmsAdults <- justTpms[, 1:9]

colSums(tpmsAdults) %>% barplot()

conditions1 <- stringr::str_extract(colnames(tpmsAdults), '[.][a-zA-Z0-9]*[.]repl')
bact <- substr(conditions1, 2, nchar(conditions1) - 5)

colnames(tpmsAdults) <- bact

head(tpmsAdults)

ttpms <- t(tpmsAdults) 
ttpmsRownames <- rownames(ttpms)
ttpmDf <- as.data.frame(ttpms)
ttpmDf$bact <- ttpmsRownames

gathered <- tidyr::gather(ttpmDf, gene, tpm, -bact)

meanTpm <- dplyr::group_by(gathered, bact, gene) %>% 
  dplyr::summarise(meanTpm = mean(tpm)) %>%
  dplyr::arrange(gene)

tpmSpread <- tidyr::spread(meanTpm, bact, meanTpm)

tpmSpread$significant <- tpmSpread$gene %in% sleuth_significant$gene

# 2629 significant genes from sleuth but only 2564 significant genes in tpmSpread

p <- ggplot2::ggplot(tpmSpread, ggplot2::aes(x = hb101, y = BIGb446, label = gene, col = significant)) +
  ggplot2::geom_point()
p
ggplot2::ggsave('output/images/scatter/tpms.pdf')
tpmSpread$gene[which(tpmSpread$hb101 == max(tpmSpread$hb101))]
### rrn-3.1
p <- p + ggplot2::scale_x_continuous(limits = c(0, 50000)) +
  ggplot2::scale_y_continuous(limits = c(0, 50000))
pPlotly <- plotly::ggplotly(p)

p <- p + ggplot2::scale_x_continuous(trans='log2') + 
  ggplot2::scale_y_continuous(trans='log2') 

library(htmlwidgets)
pPlotly <- plotly::ggplotly(p)
htmlwidgets::saveWidget(pPlotly, file.path(getwd(), paste0('output/images/scatter/tpm/', 'plotly.html')))
ggplot2::ggsave('output/images/scatter/tpm/BIGb446vsHB101.pdf')

# http://bam2rpkm.sourceforge.net/
# svn checkout http://svn.code.sf.net/p/bam2rpkm/code/trunk bam2rpkm