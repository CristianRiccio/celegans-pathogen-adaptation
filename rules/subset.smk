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
# snakemake --cores 8 --use-conda --resources downloads=1 --printshellcmds --keep-going output/genome/celegans/ref/ensembl/release95/annotation/gtf/Caenorhabditis_elegans.WBcel235Head50000.gtf

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
# snakemake --cores 8 --use-conda --resources downloads=1 --printshellcmds --keep-going output/genome/celegans/ref/ensembl/release95/annotation/gtf/Caenorhabditis_elegans.WBcel235Head{50,500,5000,800000}gene.gtf -n
