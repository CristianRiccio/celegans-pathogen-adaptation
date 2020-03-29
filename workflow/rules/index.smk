localrules: indexWormbaseGenomeWithGFF3


# rule indexBAM:
#     input:
#         bam = '{prefix}.bam'
#     output:
#         index = '{prefix}.bam.bai'
#     benchmark:
#         'output/benchmark/indexBAM/{prefix}.tsv'
#     log:
#         '{prefix}indexBAM.log'
#     conda:
#         'envs/conda/samtools1_9.yaml'
#     shell:
#         'samtools index {input.bam} > {log} 2>&1'
'''
snakemake --profile slurm output/HS678/alignment/target/celegans/ref/wormbase/release270/query/n2/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/alnSortedByCoord.bam.bai \
output/HS755/alignment/target/celegans/ref/wormbase/release270/query/nob27/embryos/hb101/repl{1,2,3}/hisat2/alnSortedByCoord.bam.bai \
output/HS755/alignment/target/celegans/ref/wormbase/release270/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/alnSortedByCoord.bam.bai \
-n
'''

rule indexHisat2bamSortedByCoord:
    input:
        bam = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bacterium}/repl{replicate}/hisat2/alnSortedByCoord.bam'
    output:
        index = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bacterium}/repl{replicate}/hisat2/alnSortedByCoord.bam.bai'
    benchmark:
        'output/benchmark/indexHisat2bamSortedByCoord/HS{HSrunNumber}/{database}/release{databaseRelease}/{wormStrain}/{stage}/{bacterium}/repl{replicate}.tsv'
    log:
        'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bacterium}/repl{replicate}/hisat2/indexHisat2bamSortedByCoord.log'
    conda:
        '../envs/conda/samtools1_9.yaml'
    shell:
        'samtools index {input.bam} > {log}'
# snakemake --cores 7 output/HS755/alignment/target/celegans/ref/ensembl/release96/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/alnSortedByCoord.bam.bai -n
'''
snakemake --profile slurm output/alignment/target/celegans/ref/ensembl/release95/query/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/alnSortedByCoord.bam.bai \
-n
'''

rule hisat2indexWormbaseGenomeWithoutGTF:
    input:
        genome = 'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/seq/c_elegans.PRJNA13758.genomic.fa'
    output:
        index = expand('output/genome/celegans/ref/wormbase/release{{wormbaseRelease}}/index/hisat2/index.{number}.ht2', number = ['1', '2', '3', '4', '5', '6', '7', '8'])
    wildcard_constraints:
        wormbaseRelease = '\d+'
    benchmark:
        'output/benchmark/hisat2indexWormbaseGenomeWithoutGTF/release{wormbaseRelease}.tsv'
    log:
        'output/log/snakemake/hisat2indexWormbaseGenomeWithoutGTF/release{wormbaseRelease}.txt'
    conda:
        '../envs/conda/hisat2v2_1_0.yaml'
    shell:
        'hisat2-build -p 1 {input.genome} output/genome/celegans/ref/wormbase/release{wildcards.wormbaseRelease}/index/hisat2/index'
# snakemake --cores 7 output/genome/celegans/ref/wormbase/release269/index/hisat2/index.{1,2,3,4,5,6,7,8}.ht2 -n

rule hisat2indexEnsemblGenomeWithoutGTF:
    input:
        genome = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/seq/genome.fa'
    output:
        index = expand('output/genome/celegans/ref/ensembl/release{{ensemblRelease}}/index/hisat2/index.{number}.ht2', number = ['1', '2', '3', '4', '5', '6', '7', '8'])
    wildcard_constraints:
        ensemblRelease = '\d+'
    benchmark:
        'output/benchmark/hisat2indexEnsemblGenomeWithoutGTF/release{ensemblRelease}.tsv'
    log:
        'output/log/snakemake/hisat2indexEnsemblGenomeWithoutGTF/release{ensemblRelease}.txt'
    conda:
        '../envs/conda/hisat2v2_1_0.yaml'
    shell:
        'hisat2-build -p 1 {input.genome} output/genome/celegans/ref/ensembl/release{wildcards.ensemblRelease}/index/hisat2/index'
# snakemake --cores 2 output/genome/celegans/ref/ensembl/release95/index/hisat2/index.{1,2,3,4,5,6,7,8}.ht2 -n

rule indexEnsemblGenomeWithGTF:
    input:
        genome = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/seq/genome.fa',
        annotation = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/annotation/gtf/Caenorhabditis_elegans.WBcel235.gtf'
    output:
        index = expand('output/genome/celegans/ref/ensembl/release{{ensemblRelease}}/index/star/withGTF/{extension}', extension = ['Genome', 'SAindex', 'SA', 'chrName.txt', 'chrStart.txt', 'chrLength.txt', 'chrNameLength.txt', 'genomeParameters.txt'])
    benchmark:
        'output/benchmark/indexEnsemblGenomeWithGTF/release{ensemblRelease}.tsv'
    log:
        'output/log/snakemake/indexEnsemblGenomeWithGTF/release{ensemblRelease}.txt'
    conda:
        '../envs/conda/star2_7_0.yaml'
    shell:
        'STAR --runMode genomeGenerate --genomeDir output/genome/celegans/ref/ensembl/release{wildcards.ensemblRelease}/index/star/withGTF --genomeFastaFiles {input.genome} --sjdbGTFfile {input.annotation} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 49 >& {log}'
# STAR --runMode genomeGenerate --genomeDir output/genome/celegans/ref/ensembl/release95/index/star/withGTF --genomeFastaFiles output/genome/celegans/ref/ensembl/release95/seq/genome.fa --sjdbGTFfile output/genome/celegans/ref/ensembl/release95/annotation/gtf/Caenorhabditis_elegans.WBcel235.gtf --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 49
# snakemake --profile slurm output/genome/celegans/ref/ensembl/release95/index/star/withGTF/Genome -n

rule indexWormbaseGenomeWithoutGFF3:
    input:
        genome = 'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/seq/c_elegans.PRJNA13758.genomic.fa'
    output:
        index = expand('output/genome/celegans/ref/wormbase/release{{wormbaseRelease}}/index/withoutGFF3/{extension}', extension = ['Genome', 'SAindex', 'SA', 'chrName.txt', 'chrStart.txt', 'chrLength.txt', 'chrNameLength.txt', 'genomeParameters.txt'])
    benchmark:
        'output/benchmark/indexWormbaseGenome/release{wormbaseRelease}.tsv'
    log:
        'output/log/conda/indexWormbaseGenome/release{wormbaseRelease}.txt'
    conda:
        '../envs/conda/star2_7_0.yaml'
    shell:
        'STAR --runMode genomeGenerate --genomeDir output/genome/celegans/ref/wormbase/release{wildcards.wormbaseRelease}/index/withoutGFF3 --genomeFastaFiles {input.genome}'
# snakemake --cores 1 output/genome/celegans/ref/wormbase/release269/index/Genome
# snakemake --profile slurm output/genome/celegans/ref/wormbase/release269/index/withoutGFF3/Genome

rule indexWormbaseGenomeWithGFF3:
    input:
        genome = 'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/seq/c_elegans.PRJNA13758.genomic.fa',
        gff3 = 'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/annotation/gff3/c_elegans.PRJNA13758.annotations.gff3'
    output:
        index = expand('output/genome/celegans/ref/wormbase/release{{wormbaseRelease}}/index/star/withGFF3/{extension}', extension = ['Genome', 'SAindex', 'SA', 'chrName.txt', 'chrStart.txt', 'chrLength.txt', 'chrNameLength.txt', 'genomeParameters.txt'])
    benchmark:
        'output/benchmark/indexWormbaseGenomeWithAnnotation/release{wormbaseRelease}.tsv'
    log:
        'output/log/snakemake/indexWormbaseGenomeWithAnnotation/release{wormbaseRelease}.txt'
    conda:
        '../envs/conda/star2_7_0b.yaml'
    shell:
        'STAR --runMode genomeGenerate --genomeDir output/genome/celegans/ref/wormbase/release{wildcards.wormbaseRelease}/index/star/withGFF3 --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gff3} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 49 >& {log}'
# snakemake --cores 1 output/genome/celegans/ref/wormbase/release269/index/star/withGFF3/Genome
'''
snakemake --profile slurm output/genome/celegans/ref/wormbase/release269/index/withGFF3/Genome
'''

rule indexTranscriptome:
    input:
        transcriptome = 'output/transcriptome/celegans/ref/{database}/release{databaseRelease}/seq/transcriptome.fa'
    output:
        index = 'output/transcriptome/celegans/ref/{database}/release{databaseRelease}/index/kallisto/transcriptome.idx'
    wildcard_constraints:
        database = '(ensembl)|(wormbase)',
        databaseRelease = '\d+',
    benchmark:
        'output/benchmark/indexTranscriptome/{database}/release{databaseRelease}.tsv'
    log:
        'output/transcriptome/celegans/ref/{database}/release{databaseRelease}/index/kallisto/transcriptome.log'
    conda:
        '../envs/conda/kallisto0_45_1.yaml'
    shell:
        'kallisto index --index {output.index} {input.transcriptome} > {log}'
# local
# snakemake --cores 7 output/transcriptome/celegans/ref/ensembl/release{94,95,96}/index/kallisto/transcriptome.idx output/transcriptome/celegans/ref/wormbase/release{268..270}/index/kallisto/transcriptome.idx -n
