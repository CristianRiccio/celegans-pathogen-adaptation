# C. elegans transgenerationally adapts to infection by regulating cysteine synthases in progeny
This repository contains all the code that was used by me (Cristian Riccio) to generate the plots and tables in Burton, Riccio et al. 2019 currently on Bioarxiv here:

Burton, Riccio et al., "Cysteine synthases CYSL-1 and CYSL-2 mediate C. elegans heritable adaptation to P. vranovensis infection", Nature Communications, 2020

https://www.nature.com/articles/s41467-020-15555-8

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7142082/

## Methods
Cutadapt version 1.18 was used to remove adapter sequences (AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC for HS678 and AGATCGGAAGAGCACACGTCTGAACTCCAGTCA (forward) and AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT (reverse) for HS755. Cutadapt was also used to trim the 3' ends of reads when the phred quality value was below 20. Reads shorter than 40 bp were discarded. The genome sequence in FASTA format and annotation file in GTF format of *C. elegans* were downloaded from Ensembl release 96. The genome sequence was indexed with the annotation with hisat2 version 2.1.0. hisat2 was used for aligning reads to the reference genome and the maximum number of alignments to be reported per read was set to 5000. Featurecounts from the conda subread package version 1.6.3 was used to count the number of reads per gene. The Ensembl release 96 annotation file in GTF format for C. elegans was used. Fragments were counted when they overlapped an exon by at least 1 nucleotide and fragments are reported at the gene level. The option for the stranded protocol was turned off. Only read pairs that had both ends aligned were counted. Given our average fragment length of 300 bp, a distance of between 50 and 600 nucleotides was tolerated for read pairs. Read pairs that had their two ends mapping to different chromosomes or mapping to the same chromosome but on different strands were not counted. Multi-mapping reads were not counted. The raw counts table was imported into R 3.5.1 for differential expression analysis with DESeq2 version 1.22.1. Normalisation was carried out within each contrast. The PCA plots were produced with DESeq2 function plotPCA() after variance stabilising transformation of the data. A snakemake workflow (CITE snakemake) was created for the RNA-seq analysis and can it be found at https://github.com/cristianriccio/celegans-pathogen-adaptation.

## Data availability
Raw sequencing data are available at the European Nucleotide Archive under accession PRJEB32993.
The URL is https://www.ebi.ac.uk/ena/browser/view/PRJEB32993.