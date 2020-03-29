Rscript code/02move.R

### Get the feature counts
snakemake --use-conda --cores 6 --printshellcmds output/counts/ensembl/featurecounts/counts.txt output/counts/ensembl/featurecounts/{adults,embryos}/{hb101,BIGb446,BIGb468}/{repl1,repl2,repl3}/counts.txt

### snakemake on the cluster
snakemake --jobs 200 --use-conda output/counts/wormbase/kallisto/{adults,embryos}/{hb101,BIGb446,BIGb468}/{repl1,repl2,repl3}/abundance.tsv --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}"

snakemake --jobs 200 output/counts/ensembl/featurecounts/{adults,embryos}/{hb101,BIGb446,BIGb468}/{repl1,repl2,repl3}/counts.txt --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}"

snakemake --jobs 200 --use-conda output/images/heatmap/adultCorrelations.pdf --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}"