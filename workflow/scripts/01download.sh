### Download data from iRODS to lustre

# 5569STDY7708519 BIGb248
# imeta qu -z seq -d sample = '5569STDY7708519'
iget -K -f /seq/pacbio/r54316_20190125_142157/2_B01/lima_output.bc1008_BAK8A--bc1008_BAK8A.bam output/reads/dna/pacbio/BIGb248/a.subreads.bam
# 5569STDY7708520 BIGb306
# 5569STDY7708521 BIGb446
# 5569STDY7708522 BIGb468
iget -K -f /seq/pacbio/r54316_20190125_142157/2_B01/m54316_190125_225603.subreads.bam output/reads/dna/pacbio/all/a.subreads.bam
