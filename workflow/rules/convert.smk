rule sam2bamHisat2:
    input:
        sam = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bact}/repl{repl}/hisat2/aln.sam'
    output:
        bam = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bact}/repl{repl}/hisat2/aln.bam'
    wildcard_constraints:
        database = '(ensembl)|(wormbase)',
        databaseRelease = '\d+',
    benchmark:
        'output/benchmark/sam2bamHisat2sam2/HS{HSrunNumber}/{database}/release{databaseRelease}/{wormStrain}/{stage}/{bact}/repl{repl}.tsv'
    log:
        'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bact}/repl{repl}/hisat2/sam2bamHisat2.log'
    conda:
        '../envs/conda/samtools1_9.yaml'
    shell:
        'samtools view -bS {input.sam} > {output.bam} 2> {log}'
# snakemake --cores 7 output/alignment/target/celegans/ref/wormbase/release95/query/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/aln.bam
# snakemake --profile slurm output/alignment/target/celegans/ref/wormbase/release95/query/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/aln.bam -n
'''
SLURM Ensembl + Wormbase and HS678 + HS755
snakemake --profile slurm output/HS678/alignment/target/celegans/ref/wormbase/release{269,270}/query/n2/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/aln.bam \
output/HS755/alignment/target/celegans/ref/wormbase/release270/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/aln.bam \
output/HS755/alignment/target/celegans/ref/wormbase/release{269,270}/query/nob27/embryos/hb101/repl{1,2,3}/hisat2/aln.bam \
output/HS678/alignment/target/celegans/ref/ensembl/release{95,96}/query/n2/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/aln.bam \
output/HS755/alignment/target/celegans/ref/ensembl/release{95,96}/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/aln.bam \
output/HS755/alignment/target/celegans/ref/ensembl/release{95,96}/query/nob27/embryos/hb101/repl{1,2,3}/hisat2/aln.bam \
-n
'''
