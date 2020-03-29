localrules: uncompressCelegansGenomeEnsembl, uncompressCelegansGFF3Ensembl, uncompressCelegansGTFensembl, uncompressCelegansGenomeWormbase, uncompressGFF3wormbase, uncompressGFF2, uncompressGTF

# C. elegans Ensembl

rule uncompressCelegansGenome:
    input:
        compressed = 'output/genome/celegans/ref/{database}/release{databaseRelease}/seq/genome.fa.gz'
    output:
        uncompressed = 'output/genome/celegans/ref/{database}/release{databaseRelease}/seq/genome.fa'
    wildcard_constraints:
        wormbaseRelease = '\d+'
    benchmark:
        'output/benchmark/uncompressCelegansGenome/{database}/release{databaseRelease}.tsv'
    log:
        'output/genome/celegans/ref/{database}/release{databaseRelease}/seq/uncompressCelegansGenome.log'
    shell:
        'gzip {input.compressed} --stdout --decompress --force > {output.uncompressed}'
# snakemake --cores 7 output/genome/celegans/ref/ensembl/release{90..96}/seq/genome.fa -n

rule uncompressCelegansTranscriptome:
    input:
        compressed = 'output/transcriptome/celegans/ref/{database}/release{databaseRelease}/seq/transcriptome.fa.gz'
    output:
        uncompressed = 'output/transcriptome/celegans/ref/{database}/release{databaseRelease}/seq/transcriptome.fa'
    wildcard_constraints:
        database = '(ensembl)|(wormbase)',
        wormbaseRelease = '\d+',
    benchmark:
        'output/benchmark/uncompressCelegansTranscriptome/{database}/release{databaseRelease}.tsv'
    log:
        'output/transcriptome/celegans/ref/{database}/release{databaseRelease}/seq/uncompressCelegansTranscriptome.stderr'
    shell:
        'gzip {input.compressed} --stdout --decompress --force > {output.uncompressed} 2> {log}'
# snakemake --cores 7 output/transcriptome/celegans/ref/ensembl/release{94,95,96}/seq/transcriptome.fa -n
# snakemake --cores 7 output/transcriptome/celegans/ref/wormbase/release{269,270}/seq/transcriptome.fa -n

rule uncompressCelegansGTFensembl:
    input:
        compressed = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/annotation/gtf/Caenorhabditis_elegans.WBcel235.gtf.gz'
    output:
        uncompressed = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/annotation/gtf/Caenorhabditis_elegans.WBcel235.gtf'
    wildcard_constraints:
        ensemblRelease = '\d+'
    benchmark:
        'output/benchmark/uncompressCelegansGTFensembl/release{ensemblRelease}.tsv'
    log:
        'output/log/snakemake/uncompressCelegansGTFensembl/release{ensemblRelease}.txt'
    shell:
        'gzip {input.compressed} --stdout --decompress --force > {output.uncompressed}'
# snakemake --cores 7 ooutput/genome/celegans/ref/ensembl/release{94,95}/annotation/gtf/Caenorhabditis_elegans.WBcel235.gtf

rule uncompressCelegansGFF3Ensembl:
    input:
        compressed = rules.downloadCelegansGFF3Ensembl.output.compressed
    output:
        uncompressed = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/annotation/gff3/Caenorhabditis_elegans.WBcel235.gff3'
    wildcard_constraints:
        wormbaseRelease = '\d+'
    benchmark:
        'output/benchmark/uncompressCelegansGFF3Ensembl/release{ensemblRelease}.tsv'
    log:
        'output/log/snakemake/uncompressCelegansGFF3Ensembl/release{ensemblRelease}.txt'
    shell:
        'gzip {input.compressed} --stdout --decompress --force > {output.uncompressed}'
# snakemake --cores 7 output/genome/celegans/ref/ensembl/release{92,93,94,95}/annotation/gff3/Caenorhabditis_elegans.WBcel235.gff3 -n

# C. elegans Wormbase

rule uncompressCelegansGenomeWormbase:
    input:
        compressed = 'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/seq/c_elegans.PRJNA13758.genomic.fa.gz'
    output:
        uncompressed = 'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/seq/c_elegans.PRJNA13758.genomic.fa'
    wildcard_constraints:
        wormbaseRelease = '\d+'
    benchmark:
        'output/benchmark/uncompressCelegansGenomeWormbase/release{wormbaseRelease}.tsv'
    log:
        'output/log/snakemake/uncompressCelegansGenomeWormbase/release{wormbaseRelease}.txt'
    shell:
        'gzip {input.compressed} --stdout --decompress --force > {output.uncompressed}'
# snakemake --cores 7 output/genome/celegans/ref/wormbase/release{260,261,262,263,264,265,266,267,268,269,270}/seq/c_elegans.PRJNA13758.genomic.fa -n

rule uncompressGFF3wormbase:
    input:
        compressed = 'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/annotation/gff3/c_elegans.PRJNA13758.annotations.gff3.gz'
    output:
        uncompressed = 'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/annotation/gff3/c_elegans.PRJNA13758.annotations.gff3'
    benchmark:
        'output/benchmark/uncompressGFF3wormbase/release{wormbaseRelease}.tsv'
    log:
        'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/annotation/gff3/uncompressGFF3wormbase.log'
    shell:
        'gzip {input.compressed} --stdout --decompress --force > {output.uncompressed}'
# snakemake --cores 7 output/genome/celegans/ref/wormbase/release{260,269,270}/annotation/gff3/c_elegans.PRJNA13758.annotations.gff3 -n

rule uncompressGFF2:
    input:
        compressed = 'output/genome/wormbase/release{wormbaseRelease}/annotation/c_elegans.PRJNA13758.annotations.gff2.gz'
    output:
        uncompressed = 'output/genome/wormbase/release{wormbaseRelease}/annotation/c_elegans.PRJNA13758.annotations.gff2'
    benchmark:
        'output/benchmark/uncompressGFF2/release{wormbaseRelease}.tsv'
    log:
        'output/log/snakemake/uncompressGFF2/release{wormbaseRelease}.txt'
    shell:
        'gzip {input.compressed} --stdout --decompress --force > {output.uncompressed}'

rule uncompressGTF:
    input:
        compressed = 'output/genome/ensembl/release{ensemblRelease}/annotation/Caenorhabditis_elegans.WBcel235.gtf.gz'
    output:
        uncompressed = 'output/genome/ensembl/release{ensemblRelease}/annotation/Caenorhabditis_elegans.WBcel235.gtf'
    benchmark:
        'output/benchmark/uncompressGTF/release{ensemblRelease}.tsv'
    log:
        'output/log/snakemake/uncompressGTF/release{ensemblRelease}.txt'
    shell:
        'gzip {input.compressed} --stdout --decompress --force > {output.uncompressed}'
# snakemake output/genome/ensembl/release95/annotation/Caenorhabditis_elegans.WBcel235.gtf -n
