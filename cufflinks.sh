#!/bin/bash
#SBATCH -n 12                   
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -p IACT
#SBATCH -o cufflinks.o
#SBATCH -e cufflinks.e
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=cr517@cam.ac.uk 
source activate cufflinks2_2_1
cufflinks -p 12 --library-type fr-firststrand --GTF output/genome/celegans/ref/ensembl/release95/annotation/gtf/Caenorhabditis_elegans.WBcel235.gtf output/alignment/target/celegans/ref/ensembl/release95/query/adults/hb101/repl1/hisat2/alnSortedByCoord.bam
