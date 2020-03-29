rule createEnsemblID_GOmap:
    input:
        script = 'workflow/scripts/createEnsemblID_GOmap.R',
    output:
        geneTable = 'output/GO/{biomartDataset}/geneTable.csv',
        ensembl_gene_id_go_id_table = 'output/GO/{biomartDataset}/ensembl_gene_id_go_id_table.csv',
        idGoMap = 'output/GO/{biomartDataset}/idGoMap.rda',
    benchmark:
        'output/benchmark/createEnsemblID_GOmap/{biomartDataset}.tsv'
    log:
        'output/GO/{biomartDataset}/createEnsemblID_GOmap.log'
    conda:
        '../envs/conda/r=3.6_bioconductor-biomart=2.42.0_r-tidyr=1.0.2.yaml'
    script:
        '../scripts/createEnsemblID_GOmap.R'
# snakemake output/GO/{celegans,mmusculus}_gene_ensembl/ensembl_gene_id_go_id_table.csv -n
