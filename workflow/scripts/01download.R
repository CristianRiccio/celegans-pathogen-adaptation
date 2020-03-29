### Download data from iRODS, rename it and transform it into fq.gz
### PacBio

snakemake --jobs 4 output/reads/dna/pacbio/{BIGb248,BIGb306,BIGb446,BIGb468}/a.subreads.bam --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" 

# 5569STDY7708519 BIGb248
imeta qu -z seq -d sample = '5569STDY7708519'
imeta ls -d /seq/pacbio/r54316_20190125_142157/2_B01/lima_output.bc1008_BAK8A--bc1008_BAK8A.bam
# library name: DN546911R-A1

# 5569STDY7708520 BIGb306
imeta qu -z seq -d sample = '5569STDY7708520'
imeta ls -d /seq/pacbio/r54316_20190125_142157/2_B01/lima_output.bc1009_BAK8A--bc1009_BAK8A.bam
# DN546911R-B1

# 5569STDY7708521 BIGb446
imeta qu -z seq -d sample = '5569STDY7708521'
imeta ls -d /seq/pacbio/r54316_20190125_142157/2_B01/lima_output.bc1010_BAK8A--bc1010_BAK8A.bam
# DN546911R-C1

# 5569STDY7708522 BIGb468
imeta qu -z seq -d sample = '5569STDY7708522'
imeta ls -d /seq/pacbio/r54316_20190125_142157/2_B01/lima_output.bc1012_BAK8A--bc1012_BAK8A.bam
# DN546911R-D1