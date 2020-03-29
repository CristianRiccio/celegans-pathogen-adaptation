localrules: compressFastaPacbio

rule compressFastaPacbio:
    input:
        fasta = 'output/reads/dna/pacbio/{pathogenORall}/a.subreads.fq'
    output:
        fastaCompressed = 'output/reads/dna/pacbio/{pathogenORall}/a.subreads.fa.gz'
    log:
        'output/log/conda/compressFastqPacbio{pathogenORall}.txt'
    shell:
        'gzip {input.fasta}'
# snakemake --profile slurm output/reads/dna/pacbio/all/a.subreads.fa.gz --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}"
