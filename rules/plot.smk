rule readsMappedPaired:
    input:
        readsMappedPaired = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/readsMappedPaired.csv',
    output:
        plot = 'output/HS{HSrunNumber}/images/bar/celegans/ref/{database}/release{databaseRelease}/nbMappedPairedReads.pdf',
    benchmark:
        'output/benchmark/readsMappedPaired/HS{HSrunNumber}/{database}/{databaseRelease}.tsv'
    conda:
        '../envs/conda/r-ggplot2.yaml'
    script:
        '../scripts/plot/readsMappedPaired.R'
# local
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS678/images/bar/celegans/ref/ensembl/release96/nbMappedPairedReads.pdf output/HS755/images/bar/celegans/ref/ensembl/release96/nbMappedPairedReads.pdf -n
# SLURM
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS678/images/bar/celegans/ref/ensembl/release96/nbMappedPairedReads.pdf output/HS755/images/bar/celegans/ref/ensembl/release96/nbMappedPairedReads.pdf --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n

rule plotCounts:
    input:
        conditions = 'samples/conditions.tsv',
        counts = 'output/HS{HSrunNumber}/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsRenamed.txt',
    output:
        summary = 'output/HS{HSrunNumber}/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsSummary.txt',
        countsTable = 'output/HS{HSrunNumber}/tables/celegans/totalCounts.csv',
        totalCounts = expand('output/HS{{HSrunNumber}}/images/bar/celegans/totalCounts.{format}', format = ['pdf', 'png', 'html']),
        boxplotCounts = expand('output/HS{{HSrunNumber}}/images/boxplot/counts.{format}', format = ['pdf', 'png', 'html']),
        boxplotCountsLog10 = expand('output/HS{{HSrunNumber}}/images/boxplot/countsLog10.{format}', format = ['pdf', 'png', 'html']),
        histogramCounts = expand('output/HS{{HSrunNumber}}/images/histogram/counts.{format}', format = ['pdf', 'png', 'html']),
        histogramCountsLog10 = expand('output/HS{{HSrunNumber}}/images/histogram/countsLog10.{format}', format = ['pdf', 'png', 'html']),
    benchmark:
        'output/benchmark/plotCounts/HS{HSrunNumber}.tsv'
    conda:
        '../envs/conda/r-ggplot2_r-plotly_r-stringr_r-tidyr.yaml'
    script:
        '../scripts/plot/plotCounts.R'
# local
'''
snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason \
output/HS{678,755}/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsSummary.txt \
output/HS{678,755}/tables/celegans/totalCounts.csv \
output/HS{678,755}/images/bar/celegans/totalCounts.{pdf,png,html} \
output/HS{678,755}/images/boxplot/counts.{pdf,png,html} \
output/HS{678,755}/images/boxplot/countsLog10.{pdf,png,html} \
output/HS{678,755}/images/histogram/counts.{pdf,png,html} \
output/HS{678,755}/images/histogram/countsLog10.{pdf,png,html} \
-n
'''
# SLURM
'''
snakemake --jobs 2000 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason \
output/HS{678,755}/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsSummary.txt \
output/HS{678,755}/tables/celegans/totalCounts.csv \
output/HS{678,755}/images/bar/celegans/totalCounts.{pdf,png,html} \
output/HS{678,755}/images/boxplot/counts.{pdf,png,html} \
output/HS{678,755}/images/boxplot/countsLog10.{pdf,png,html} \
output/HS{678,755}/images/histogram/counts.{pdf,png,html} \
output/HS{678,755}/images/histogram/countsLog10.{pdf,png,html} \
--cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} \
-t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} \
--error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} \
--mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n
'''
