rule validate_fastq:
    input:
        '{filename}.fq.gz',
    output:
        '{filename}_fqtools_validate.txt',
    conda:
        '../envs/conda/fqtools=2.0.yaml'
    shell:
        'fqtools validate {input} > {output}'
'''
snakemake --cores 2 --use-conda --keep-going output/HS678/reads/rna/illumina/celegans/merged/n2/{adults,embryos}/{hb101,BIGb446,BIGb468}/repl{1,2,3}/R{1,2}/L1L2_fqtools_validate.txt \
output/HS755/reads/rna/illumina/celegans/merged/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/R{1,2}/L1L2_fqtools_validate.txt \
output/HS755/reads/rna/illumina/celegans/merged/nob27/embryos/hb101/repl{1,2,3}/R{1,2}/L1L2_fqtools_validate.txt \
-n
'''
