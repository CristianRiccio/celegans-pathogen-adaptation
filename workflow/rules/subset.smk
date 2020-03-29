rule headGTF:
    input:
        'output/genome/celegans/ref/ensembl/release95/annotation/gtf/Caenorhabditis_elegans.WBcel235.gtf'
    output:
        'output/genome/celegans/ref/ensembl/release95/annotation/gtf/Caenorhabditis_elegans.WBcel235Head{subsetSize}.gtf'
    wildcard_constraints:
        subsetSize = '\d*',
    benchmark:
        'output/benchmark/subset/subsetSize{subsetSize}.tsv'
    log:
        'output/log/snakemake/subset/{subsetSize}.txt'
    shell:
        'head -n {wildcards.subsetSize} {input} > {output} 2> {log}'
# snakemake --cores 8 output/genome/celegans/ref/ensembl/release95/annotation/gtf/Caenorhabditis_elegans.WBcel235Head50000.gtf -n

rule subsetGTFbyFeature:
    input:
        gtf = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/annotation/gtf/Caenorhabditis_elegans.WBcel235Head{subsetSize}.gtf'
    output:
        gtfFeature = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/annotation/gtf/Caenorhabditis_elegans.WBcel235Head{subsetSize}gene.gtf'
    wildcard_constraints:
        ensemblRelease = '\d+',
        subsetSize = '\d+',
    benchmark:
        'output/benchmark/subsetGTFbyFeature/release{ensemblRelease}subsetSize{subsetSize}.tsv'
    log:
        'output/log/snakemake/subsetGTFbyFeature/release{ensemblRelease}subsetSize{subsetSize}.txt'
    conda:
        '../envs/conda/python3_7_3.yaml'
    script:
        '../scripts/subsetGTFbyFeature.py'
# snakemake --cores 8 output/genome/celegans/ref/ensembl/release95/annotation/gtf/Caenorhabditis_elegans.WBcel235Head{50,500,5000,800000}gene.gtf -n

rule onlyKeepCounts:
    input:
        script = 'workflow/scripts/subset/onlyKeepCounts.R',
        counts = '{prefix}Counts.txt',
    output:
        countsOnly = '{prefix}CountsOnly.csv',
    benchmark:
        'output/benchmark/onlyKeepCounts/{prefix}.tsv'
    log:
        '{prefix}_onlyKeepCounts.log'
    conda:
        '../../envs/conda/r=3.6.yaml'
    priority: 100
    threads: 1
    script:
        '../scripts/subset/onlyKeepCounts.R'
# snakemake --profile slurm output/HS678HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsOnly.csv -n
