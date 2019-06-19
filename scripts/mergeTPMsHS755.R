# Merge pseudocounts from the 755 run into one table

# Put the folder in which the tpms are stored into a variable
folder <- 'output/HS755/counts/ensembl/release96/kallisto/'

# Store the 9 sample names for N2 into a variable
samplesN2 <- paste0('n2_', rep(c('hb101', 'pa14', 'plum'), each = 3), '_repl', rep(c(1, 2, 3), 3))
# Store the 3 sample names for NOB27 into a variable
samplesNOB27 <- paste0('nob27_hb101_repl', c(1, 2, 3))
# Store the 9 paths to tpms for N2 into a variable
paths2abundancesN2 <- paste0(folder, 'n2', '/embryos/', rep(c('hb101', 'pa14', 'plum'), each = 3), '/repl', rep(c(1, 2, 3), 3), '/abundance.tsv')
# Store the 3 paths to tpms for NOB27 into a variable
paths2abundancesNOB27 <- paste0(folder, 'nob27/embryos/hb101/repl', c(1, 2, 3), '/abundance.tsv')
# Concatenate sample names
samples <- c(samplesN2, samplesNOB27)
# Concatenate the paths to tpms
paths2abundances <- c(paths2abundancesN2, paths2abundancesNOB27)

abundances <- read.delim(paths2abundances[1])
colnames(abundances)[5] <- samples[1]

for (i in 2:length(samples)) {
  abundancesi <- read.delim(paths2abundances[i])[, c(1, 5)]
  colnames(abundancesi)[2] <- samples[i]
  abundances <- merge(abundances, abundancesi, by.x = 'target_id', by.y = 'target_id')
}

write.csv(abundances, paste0(folder, 'tpms.csv'), 
          row.names = FALSE, quote = FALSE)

