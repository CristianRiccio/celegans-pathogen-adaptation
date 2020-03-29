localrules: downloadCelegansGenomeEnsembl, downloadCelegansTranscriptomeEnsembl
localrules: downloadCelegansGTFensembl, downloadCelegansGenomeWormbase
localrules: downloadCelegansTranscriptomeWormbase, downloadCelegansGFF3wormbase
localrules: downloadCelegansGFF2wormbase,

# C. elegans Ensembl

rule downloadCelegansGenomeEnsembl:
    output:
        compressed = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/seq/genome.fa.gz'
    benchmark:
        'output/benchmark/downloadEnsemblGenome/release{ensemblRelease}.tsv'
    log:
        'output/genome/celegans/ref/ensembl/release{ensemblRelease}/seq/downloadCelegansGenomeEnsembl.log'
    resources:
        downloads = 1
    conda:
        '../envs/conda/wget.yaml'
    shell:
        'wget ftp://ftp.ensembl.org/pub/release-{wildcards.ensemblRelease}/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz -O {output.compressed} > {log} 2>&1'
# snakemake --cores 7 output/genome/celegans/ref/ensembl/release{90..96}/seq/genome.fa.gz -n

rule downloadCelegansTranscriptomeEnsembl:
    output:
        compressed = 'output/transcriptome/celegans/ref/ensembl/release{ensemblRelease}/seq/transcriptome.fa.gz'
    benchmark:
        'output/benchmark/downloadCelegansTranscriptomeEnsembl/release{ensemblRelease}.tsv'
    log:
        'output/transcriptome/celegans/ensembl/release{ensemblRelease}/seq/downloadCelegansTranscriptomeEnsembl.log'
    resources:
        downloads = 1
    conda:
        '../envs/conda/wget.yaml'
    shell:
        'wget ftp://ftp.ensembl.org/pub/release-{wildcards.ensemblRelease}/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz -O {output.compressed}'
# snakemake --cores 7 output/transcriptome/celegans/ref/ensembl/release{94,95,96}/seq/transcriptome.fa.gz -n

rule downloadCelegansGTFensembl:
    output:
        compressed = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/annotation/gtf/Caenorhabditis_elegans.WBcel235.gtf.gz'
    benchmark:
        'output/benchmark/downloadCelegansGTFensembl/release{ensemblRelease}.tsv'
    log:
        'output/log/snakemake/downloadCelegansGTFensembl/release{ensemblRelease}.txt'
    resources:
        downloads = 1
    conda:
        '../envs/conda/wget.yaml'
    shell:
        'wget ftp://ftp.ensembl.org/pub/release-{wildcards.ensemblRelease}/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.{wildcards.ensemblRelease}.gtf.gz -O {output.compressed}'
# snakemake --cores 7 output/genome/celegans/ref/ensembl/release{94,95,96}/annotation/gtf/Caenorhabditis_elegans.WBcel235.gtf.gz -n

rule downloadCelegansGFF3Ensembl:
    output:
        compressed = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/annotation/gff3/Caenorhabditis_elegans.WBcel235.gff3.gz'
    benchmark:
        'output/benchmark/downloadCelegansGFF3Ensembl/release{ensemblRelease}.tsv'
    log:
        'output/log/snakemake/downloadCelegansGFF3Ensembl/release{ensemblRelease}.txt'
    conda:
        '../envs/conda/wget.yaml'
    resources:
        downloads = 1
    shell:
        'wget ftp://ftp.ensembl.org/pub/release-{wildcards.ensemblRelease}/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.{wildcards.ensemblRelease}.gff3.gz -O {output.compressed} > {log} 2>&1'
# snakemake --cores 7 output/genome/celegans/ref/ensembl/release{92,93,94,95}/annotation/gff3/Caenorhabditis_elegans.WBcel235.gff3.gz -n

# C. elegans Wormbase

rule downloadCelegansGenomeWormbase:
    output:
        compressed = 'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/seq/c_elegans.PRJNA13758.genomic.fa.gz'
    benchmark:
        'output/benchmark/downloadCelegansGenomeWormbase/release{wormbaseRelease}.tsv'
    log:
        'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/seq/downloadCelegansGenomeWormbase.log'
    resources:
        downloads = 1
    conda:
        '../envs/conda/wget.yaml'
    shell:
        'wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/genomic/c_elegans.PRJNA13758.WS{wildcards.wormbaseRelease}.genomic.fa.gz -O {output.compressed} > {log} 2>&1'
# snakemake --cores 7 output/genome/celegans/ref/wormbase/release{260,261,262,263,264,265,266,267,268,269,270}/seq/c_elegans.PRJNA13758.genomic.fa.gz -n

rule downloadCelegansTranscriptomeWormbase:
    output:
        compressed = 'output/transcriptome/celegans/ref/wormbase/release{wormbaseRelease}/seq/transcriptome.fa.gz'
    benchmark:
        'output/benchmark/downloadCelegansTranscriptomeWormbase/release{wormbaseRelease}.tsv'
    log:
        'output/transcriptome/celegans/ref/wormbase/release{wormbaseRelease}/seq/downloadCelegansTranscriptomeWormbase.log'
    conda:
        '../envs/conda/wget.yaml'
    shell:
        'wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.WS{wildcards.wormbaseRelease}.mRNA_transcripts.fa.gz -O {output.compressed}'
# local
# snakemake --cores 7 output/transcriptome/celegans/ref/wormbase/release{269,270}/seq/transcriptome.fa.gz -n

rule downloadCelegansGFF3wormbase:
    output:
        compressed = 'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/annotation/gff3/c_elegans.PRJNA13758.annotations.gff3.gz'
    benchmark:
        'output/benchmark/downloadCelegansGFF3wormbase/release{wormbaseRelease}.tsv'
    log:
        'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/annotation/gff3/downloadCelegansGFF3wormbase.log'
    resources:
        downloads = 1
    conda:
        '../envs/conda/wget.yaml'
    shell:
        'wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/gff/c_elegans.PRJNA13758.WS{wildcards.wormbaseRelease}.annotations.gff3.gz -O {output.compressed}'
# snakemake --cores 7 output/genome/celegans/ref/wormbase/release{260,269,270}/annotation/gff3/c_elegans.PRJNA13758.annotations.gff3.gz -n

rule downloadCelegansGFF2wormbase:
    output:
        compressed = 'output/genome/celegans/wormbase/release{wormbaseRelease}/annotation/c_elegans.PRJNA13758.annotations.gff2.gz'
    benchmark:
        'output/benchmark/downloadCelegansGFF2wormbase/release{wormbaseRelease}.tsv'
    log:
        'output/log/conda/downloadCelegansGFF2wormbase/release{wormbaseRelease}.txt'
    resources:
        downloads = 1
    conda:
        '../envs/conda/wget.yaml'
    shell:
        'wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/gff/c_elegans.PRJNA13758.WS{wildcards.wormbaseRelease}.annotations.gff2.gz -O {output.compressed}'
# snakemake --cores 7 output/genome/celegans/wormbase/release{260,269}/annotation/c_elegans.PRJNA13758.annotations.gff2.gz -n
