localrules: moveAndRenameData, mergeLanes, nbMappedPairedReads

import pandas
samples = pandas.read_csv('samples/datasets.tsv', sep = '\t')

# Transfer data between the Gurdon cluster and my laptop
# rsync -avz cr517@cb-head2.gurdon.private.cam.ac.uk:/mnt/home2/miska/cr517/projects/celegans-pathogen-adaptation/code/output/HS678HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt /Users/cr517/Documents/phd/projects/celegans-pathogen-adaptation/code/output/HS678HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt -n

# Transfer data between my laptop and the Gurdon cluster
# rsync -avz /Users/cr517/Documents/phd/projects/celegans-pathogen-adaptation/code/{config.yaml,envs/createCondaEnvs.sh,input,rules,samples,scripts,slurm.json,Snakefile} cr517@cb-head2.gurdon.private.cam.ac.uk:/mnt/home2/miska/cr517/projects/celegans-pathogen-adaptation/code/ -n
# rsync -avz /Users/cr517/Documents/phd/projects/celegans-pathogen-adaptation/code/envs/createCondaEnvs.sh cr517@cb-head2.gurdon.private.cam.ac.uk:/mnt/home2/miska/cr517/projects/celegans-pathogen-adaptation/code/envs/ -n

# mkdir -p output/images/graph/rulegraph/
# Create graphs to visualise the pipeline
# snakemake --rulegraph output/counts/featurecounts/hisat2/celegans/ensembl/release95/minoverlap1strandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600DiscardMisalignedPairsNoMultimappersFPM.txt | dot -Tpdf > output/images/graph/rulegraph/counts.pdf
# snakemake --rulegraph | dot -Tpdf > output/images/graph/rulegraph/target.pdf
# snakemake --dag | dot -Tpdf > output/images/graph/dag/target.pdf

configfile: "config.yaml"

# On 25th May 2018, the latest release of Ensembl is 96
# Define variables
ensemblRelease = '96'

coresPerNode = 24
stages = ['adults', 'embryos']
bacteria = ['hb101', 'BIGb446', 'BIGb468', 'pa14', 'plum']
replicates = [1, 2, 3]

# Experiment 2
n2bacteria = ['hb101', 'pa14', 'plum']
nob27bacterium = ['hb101']

# Create regexes for stages, bacteria and replicates
stagesRE = '(adults)|(embryos)'
bacteriaRE = '(hb101)|(BIGb446)|(BIGb468)|(pa14)|(plum)'
replicatesRE = '(repl1)|(repl2)|(repl3)'

# import python modules
import os
newpath = r'output/log/job/error'
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = r'output/log/job/output'
if not os.path.exists(newpath):
    os.makedirs(newpath)

include: 'rules/download.smk'
include: 'rules/compressUncompress.smk'
include: 'rules/index.smk'
include: 'rules/align.smk'
include: 'rules/convert.smk'
include: 'rules/subset.smk'
include: 'rules/plot.smk'

rule target:
    input:
        'output/transcriptome/celegans/ref/' + str(config['database']) + '/release' + str(config['databaseRelease']) + '/index/kallisto/transcriptome.idx',
        expand('output/HS{HSrunNumber}/counts/featurecounts/hisat2/celegans/' + str(config['database']) + '/release' + str(config['databaseRelease']) + '/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersFPKMs.txt', HSrunNumber = [678, 755]),
        expand('output/HS{HSrunNumber}/images/bar/celegans/ref/' + str(config['database']) + '/release' + str(config['databaseRelease']) + '/nbMappedPairedReads.pdf', HSrunNumber = [678, 755]),
        # expand('output/genome/wormbase/release267/index/{extension}', extension = ['Genome', 'SAindex', 'SA', 'chrName.txt', 'chrStart.txt', 'chrLength.txt', 'chrNameLength.txt', 'genomeParameters.txt']),
        # 'output/genome/wormbase/release267/annotation/c_elegans.PRJNA13758.annotations.gff3',
        # 'output/genome/wormbase/release267/annotation/c_elegans.PRJNA13758.annotations.gff2',
        # expand('output/counts/wormbase/release267/kallisto/{stage}/{bact}/{repl}/abundance.tsv', stage = stages, bact = bacteria, repl = replicates),
        # expand('output/alignment/ensembl/release95/{stage}/{bact}/{repl}/Aligned.sortedByName.out.bam', stage = stages, bact = bacteria, repl = replicates),
        # expand('output/counts/featurecounts/ensembl/release95/{stage}/{bact}/{repl}/counts.txt',  stage = stages, bact = bacteria, repl = replicates),
        # 'output/counts/featurecounts/ensembl/release' + ensemblRelease + '/countsMultimappersMultioverlappers.txt',
        # expand('output/counts/ensembl/release95/htseq/{stage}/{bact}/{repl}/counts.txt', stage = stages, bact = bacteria, repl = replicates),
        # 'output/robjects/featurecounts/ensembl/release95/dds2.Rda',
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going -n
# SLURM
# snakemake --jobs 200 --resources downloads=1 --restart-times 1 --use-conda --printshellcmds --keep-going --reason --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n

# Copy the raw data from the input folder to the output folder
# I downloaded the raw data from https://basespace.illumina.com
# I then moved the datasets out of their subfolders into
# the folder named HS678 (or HS755) using `mv input/HS755/*/* input/HS755`
# (or `mv input/HS678/*/* input/HS678`)
# I then manually deleted all the empty folders from the Finder
# Note that it is OK to move the datasets out of their subfolders
# and delete the subfolders because the subfolder names
# contain a subset of the information contained in the dataset names
# Then, transfer the raw data to SLURM

def findDataset(wildcards):
    import pandas
    HSrunNumber = int(wildcards[0])
    wormStrain = wildcards[1]
    stage = wildcards[2]
    bacterium = wildcards[3]
    replicate = int(wildcards[4])
    readOrientation = int(wildcards[5])
    lane = int(wildcards[6])
    samples = pandas.read_csv('samples/datasets.tsv', sep = '\t')
    dataset = samples.dataset[(samples.HSrunNumber == HSrunNumber) & (samples.wormStrain == wormStrain) & (samples.stage == stage) & (samples.bacterium == bacterium) & (samples.replicate == replicate) & (samples.lane == lane) & (samples.readOrientation == readOrientation)]
    # File path 'input/HS678/N2-adult-HB101-3_S3_L002_R2_001.fastq.gz  '
    # ends with whitespace.
    # This is likely unintended.
    # It can also lead to inconsistent results of the file-matching
    # approach used by Snakemake.
    # That's why I used .rstrip() to remove trailing whitespace at the end
    # of the name of dataset. This is a more robust approach than manually
    # removing whitespace from the datasets table.
    return 'input/HS' + str(HSrunNumber) + '/' + str(dataset.iloc[0]).rstrip()

rule moveAndRenameData:
    input: findDataset
    output:
        temp('output/HS{HSrunNumber}/reads/rna/illumina/celegans/raw/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R{readOrientation}/lane/L{lane}.fq.gz')
    shell:
        'cp {input} {output}'
# local
# snakemake --cores 8 --resources downloads=1 irods=10 --use-conda --printshellcmds --keep-going --reason output/HS678/reads/rna/illumina/celegans/raw/n2/{adults,embryos}/{hb101,BIGb446,BIGb468}/repl{1,2,3}/R{1,2}/lane/L{1,2}.fq.gz output/HS755/reads/rna/illumina/celegans/raw/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/R{1,2}/lane/L{1,2}.fq.gz output/HS755/reads/rna/illumina/celegans/raw/nob27/embryos/hb101/repl{1,2,3}/R{1,2}/lane/L{1,2}.fq.gz -n

# Merge fastq files together that come from the
# same sample and have reads in the same direction
# but have been produced in different lanes
# (The sequencing technician from the Gurdon Institute showed me the 2 lanes on the machine)
rule mergeLanes:
    input:
        fastq1 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/raw/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R{readOrientation}/lane/L1.fq.gz',
        fastq2 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/raw/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R{readOrientation}/lane/L2.fq.gz',
    output:
        temp('output/HS{HSrunNumber}/reads/rna/illumina/celegans/merged/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R{readOrientation}/L1L2.fq.gz')
    benchmark:
        'output/benchmark/mergeLanes/HS{HSrunNumber}/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R{readOrientation}.tsv'
    shell:
        'cat {input.fastq1} {input.fastq2} > {output}'
# local
# snakemake --cores 8 --resources downloads=1 irods=10 --use-conda --printshellcmds --keep-going --reason output/HS678/reads/rna/illumina/celegans/merged/n2/{adults,embryos}/{hb101,BIGb446,BIGb468}/repl{1,2,3}/R{1,2}/L1L2.fq.gz output/HS755/reads/rna/illumina/celegans/merged/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/R{1,2}/L1L2.fq.gz output/HS755/reads/rna/illumina/celegans/merged/nob27/embryos/hb101/repl{1,2,3}/R{1,2}/L1L2.fq.gz -n

# Remove adapter sequence, low quality bases at the end of the reads (< Q20) and reads shorter than 40 bp
rule cutadapt:
    input:
        reads1 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/merged/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R1/L1L2.fq.gz',
        reads2 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/merged/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R2/L1L2.fq.gz',
    output:
        trimmed1 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/trimmed/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R1/L1L2.fq.gz',
        trimmed2 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/trimmed/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R2/L1L2.fq.gz',
        tooShort1 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/trimmed/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R1/tooShort.fa',
        tooShort2 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/trimmed/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R2/tooShort.fa',
    wildcard_constraints:
        stage = stagesRE,
        bact = bacteriaRE,
        repl = replicatesRE
    benchmark:
        'output/benchmark/cutadapt/HS{HSrunNumber}/{wormStrain}/{stage}/{bacterium}/repl{replicate}.tsv'
    log:
        'output/HS{HSrunNumber}/reads/rna/illumina/celegans/trimmed/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R1/cutadapt.log',
    params:
        adapterRead1 = lambda wildcards: config['adapterRead1HS' + str(wildcards[0])],
        adapterRead2 = lambda wildcards: config['adapterRead2HS' + str(wildcards[0])],
    conda:
        'envs/conda/cutadapt1_18.yaml'
    shell:
        'cutadapt --adapter {params.adapterRead1} -A {params.adapterRead2} -q 20 --minimum-length 40 --too-short-output {output.tooShort1} --too-short-paired-output {output.tooShort2} --output {output.trimmed1} -p {output.trimmed2} {input.reads1} {input.reads2} > {log} 2>&1'
# local
# snakemake --cores 7 --use-conda --resources downloads=1 --printshellcmds --keep-going --reason output/HS678/reads/rna/illumina/celegans/trimmed/n2/{adults,embryos}/{hb101,BIGb446,BIGb468}/repl{1,2,3}/R{1,2}/L1L2.fq.gz output/HS755/reads/rna/illumina/celegans/trimmed/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/R{1,2}/L1L2.fq.gz output/HS755/reads/rna/illumina/celegans/trimmed/nob27/embryos/hb101/repl{1,2,3}/R{1,2}/L1L2.fq.gz -n
# SLURM
# snakemake --jobs 2000 --use-conda --resources downloads=1 --printshellcmds --keep-going --reason output/HS678/reads/rna/illumina/celegans/trimmed/n2/{adults,embryos}/{hb101,BIGb446,BIGb468}/repl{1,2,3}/R{1,2}/L1L2.fq.gz output/HS755/reads/rna/illumina/celegans/trimmed/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/R{1,2}/L1L2.fq.gz output/HS755/reads/rna/illumina/celegans/trimmed/nob27/embryos/hb101/repl{1,2,3}/R{1,2}/L1L2.fq.gz --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n

# Run FASTQC on the merged reads (before adapter, quality and length trimming)
# and on the trimmed reads (after adapter, qualilty and length trimming)
rule fastqc:
    input:
        reads = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/merged/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R{readOrientation}/L1L2.fq.gz',
    output:
        QCreport = expand('output/HS{{HSrunNumber}}/reads/rna/illumina/celegans/{{reads}}/{{wormStrain}}/{{stage}}/{{bacterium}}/repl{{replicate}}/R{{readOrientation}}/qc/L1L2_fastqc.{extension}', extension = ['html', 'zip'])
    wildcard_constraints:
        stage = stagesRE,
        bact = bacteriaRE,
        repl = replicatesRE,
        readOrientation = '1|2',
    benchmark:
        'output/benchmark/fastqc/HS{HSrunNumber}/{reads}/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R{readOrientation}.tsv'
    log:
        'output/HS{HSrunNumber}/reads/rna/illumina/celegans/{reads}/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R{readOrientation}/qc/fastqc.log',
    conda:
        'envs/conda/fastqc0_11_8.yaml'
    shell:
        'fastqc {input.reads} -o output/HS{wildcards.HSrunNumber}/reads/rna/illumina/celegans/{wildcards.reads}/{wildcards.wormStrain}/{wildcards.stage}/{wildcards.bacterium}/repl{wildcards.replicate}/R{wildcards.readOrientation}/qc/ > {log} 2>&1'

# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/reads/rna/illumina/celegans/{merged,trimmed}/{adults,embryos}/{hb101,BIGb446,BIGb468}/repl{1,2,3}/qc/R{1,2}_fastqc.{html,zip}
# HS755
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason  output/HS755/reads/rna/illumina/celegans/merged/n2/embryos/hb101/repl{1,2,3}/R{1,2}/qc/L1L2_fastqc.{html,zip} -n
# snakemake --jobs 200 --resources downloads=1 --use-conda --printshellcmds --keep-going output/reads/rna/illumina/celegans/{merged,trimmed}/{adults,embryos}/{hb101,BIGb446,BIGb468}/repl{1,2,3}/qc/R{1,2}_fastqc.{html,zip} --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}"

rule multiqcExp1:
    input:
        QCreports = expand('output/reads/rna/illumina/celegans/{reads}/{stage}/{bact}/{repl}/qc/R{readOrientation}_fastqc.{extension}', reads = ['merged', 'trimmed'], stage = stages, bact = bacteria, repl = replicates, readOrientation = [1, 2], extension = ['html', 'zip'])
    output:
        'output/reads/rna/illumina/celegans/qc/multiqc_report.html',
        expand('output/reads/rna/illumina/celegans/qc/multiqc_data/{file}', file = ['multiqc_data.json', 'multiqc_fastqc.txt', 'multiqc_general_stats.txt', 'multiqc.log', 'multiqc_sources.txt']),
    benchmark:
        'output/benchmark/multiqc/benchmark.tsv'
    log:
        output = 'output/log/snakemake/output/multiqc/output.o',
        error = 'output/log/snakemake/error/multiqc/error.e',
    conda:
        'envs/conda/multiqc1_7.yaml'
    shell:
        'multiqc --force --dirs --outdir output/reads/rna/illumina/celegans/qc output/reads/rna/illumina/celegans/'
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going output/reads/rna/illumina/celegans/qc/multiqc_report.html
# snakemake --jobs 200 --resources downloads=1 --use-conda --printshellcmds --keep-going output/reads/rna/illumina/celegans/qc/multiqc_report.html --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus}--mail-type {cluster.mailType} --mail-user {cluster.mailUser}"

rule multiqcExp2:
    input:
        QCreports = expand('output/HS755/reads/rna/illumina/celegans/{reads}/n2/embryos/{bacterium}/repl{replicate}/R{readOrientation}/qc/L1L2_fastqc.{extension}', reads = ['merged', 'trimmed'], bacterium = n2bacteria, replicate = [1, 2, 3], readOrientation = [1, 2], extension = ['html', 'zip'])
    output:
        'output/HS755/reads/rna/illumina/celegans/qc/multiqc_report.html',
        expand('output/HS755/reads/rna/illumina/celegans/qc/multiqc_data/{file}', file = ['multiqc_data.json', 'multiqc_fastqc.txt', 'multiqc_general_stats.txt', 'multiqc.log', 'multiqc_sources.txt']),
    benchmark:
        'output/benchmark/multiqc/benchmark.tsv'
    log:
        'output/HS755/reads/rna/illumina/celegans/qc/multiqcExp2.log',
    conda:
        'envs/conda/multiqc1_7.yaml'
    shell:
        'multiqc --force --dirs --outdir output/HS755/reads/rna/illumina/celegans/qc output/HS755/reads/rna/illumina/celegans/ > {log} 2>&1'
# local
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS755/reads/rna/illumina/celegans/qc/multiqc_report.html -n
# snakemake --jobs 200 --resources downloads=1 --use-conda --printshellcmds --keep-going output/reads/rna/illumina/celegans/qc/multiqc_report.html --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}"

rule sortBamByCoord:
    input:
        bam = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bact}/repl{repl}/hisat2/aln.bam'
    output:
        bamSorted = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bact}/repl{repl}/hisat2/alnSortedByCoord.bam'
    wildcard_constraints:
        database = '(ensembl)|(wormbase)',
        databaseRelease = '\d+',
    benchmark:
        'output/benchmark/sortBamByCoord/HS{HSrunNumber}/{database}/release{databaseRelease}/{wormStrain}/{stage}/{bact}/repl{repl}.tsv'
    log:
        'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bact}/repl{repl}/hisat2/sortBamByCoord.log'
    conda:
        'envs/conda/samtools1_9.yaml'
    shell:
        'samtools sort -o {output.bamSorted} --output-fmt BAM {input.bam} > {log} 2>&1'
# snakemake --cores 7 --use-conda --printshellcmds --keep-going output/alignment/target/celegans/ref/ensembl/release95/query/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/alnSortedByCoord.bam
# snakemake --jobs 7 --use-conda --resources downloads=1 --printshellcmds --keep-going output/alignment/target/celegans/ref/ensembl/release95/query/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/alnSortedByCoord.bam --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n
# SLURM Ensembl + Wormbase and HS678 + HS755
# snakemake --jobs 2000 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS678/alignment/target/celegans/ref/wormbase/release270/query/n2/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/alnSortedByCoord.bam output/HS755/alignment/target/celegans/ref/wormbase/release270/query/nob27/embryos/hb101/repl{1,2,3}/hisat2/alnSortedByCoord.bam output/HS755/alignment/target/celegans/ref/wormbase/release270/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/alnSortedByCoord.bam output/HS678/alignment/target/celegans/ref/ensembl/release{95,96}/query/n2/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/alnSortedByCoord.bam output/HS755/alignment/target/celegans/ref/ensembl/release{95,96}/query/nob27/embryos/hb101/repl{1,2,3}/hisat2/alnSortedByCoord.bam output/HS755/alignment/target/celegans/ref/ensembl/release{95,96}/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/alnSortedByCoord.bam --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n

rule alignmentStats:
    input:
        aln = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bact}/repl{replicate}/hisat2/alnSortedByCoord.bam'
    output:
        stats = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bact}/repl{replicate}/hisat2/alnSortedByCoordStats.txt'
    benchmark:
        'output/benchmark/alignmentStats/HS{HSrunNumber}/{database}/release{databaseRelease}/{wormStrain}/{stage}/{bact}/repl{replicate}.tsv'
    log:
        error = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bact}/repl{replicate}/hisat2/alignmentStats.err'
    conda:
        'envs/conda/samtools1_9.yaml'
    shell:
        'samtools stats {input.aln} > {output.stats} 2> {log.error}'
# local
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS755/alignment/target/celegans/ref/ensembl/release95/query/n2/embryos/hb101/repl1/hisat2/alnSortedByCoordStats.txt -n
# SLURM Ensembl + Wormbase and HS678 + HS755
'''
snakemake --jobs 200 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason \
output/HS678/alignment/target/celegans/ref/wormbase/release270/query/n2/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/alnSortedByCoordStats.txt \
output/HS678/alignment/target/celegans/ref/ensembl/release{95,96}/query/n2/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/alnSortedByCoordStats.txt \
output/HS755/alignment/target/celegans/ref/ensembl/release{95,96}/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/alnSortedByCoordStats.txt \
output/HS755/alignment/target/celegans/ref/ensembl/release{95,96}/query/nob27/embryos/hb101/repl{1,2,3}/hisat2/alnSortedByCoordStats.txt \
output/HS755/alignment/target/celegans/ref/wormbase/release{269,270}/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/alnSortedByCoordStats.txt \
output/HS755/alignment/target/celegans/ref/wormbase/release{269,270}/query/nob27/embryos/hb101/repl{1,2,3}/hisat2/alnSortedByCoordStats.txt \
--cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} \
--error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n
'''

rule nbMappedPairedReads:
    input:
        alnStats = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bact}/repl{replicate}/hisat2/alnSortedByCoordStats.txt'
    output:
        nbMappedPairedReads = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bact}/repl{replicate}/hisat2/alnSortedByCoordStatsReadsMappedAndPaired.txt'
    benchmark:
        'output/benchmark/nbMappedPairedReads/HS{HSrunNumber}/{database}/release{databaseRelease}/{wormStrain}/{stage}/{bact}/{replicate}.tsv'
    log:
        'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/{wormStrain}/{stage}/{bact}/repl{replicate}/hisat2/nbMappedPairedReads.stderr'
    shell:
        "grep 'reads mapped and paired' {input.alnStats} | egrep -o '[0-9]*' > {output.nbMappedPairedReads} 2> {log}"
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS678/alignment/target/celegans/ref/ensembl/release{95,96}/query/n2/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/alnSortedByCoordStatsReadsMappedAndPaired.txt output/HS755/alignment/target/celegans/ref/ensembl/release{95,96}/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/alnSortedByCoordStatsReadsMappedAndPaired.txt output/HS755/alignment/target/celegans/ref/ensembl/release{95,96}/query/nob27/embryos/hb101/repl{1,2,3}/hisat2/alnSortedByCoordStatsReadsMappedAndPaired.txt -n
# SLURM
# snakemake --jobs 200 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS678/alignment/target/celegans/ref/ensembl/release{95,96}/query/n2/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/alnSortedByCoordStatsReadsMappedAndPaired.txt output/HS755/alignment/target/celegans/ref/ensembl/release{95,96}/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/alnSortedByCoordStatsReadsMappedAndPaired.txt output/HS755/alignment/target/celegans/ref/ensembl/release{95,96}/query/nob27/embryos/hb101/repl{1,2,3}/hisat2/alnSortedByCoordStatsReadsMappedAndPaired.txt --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n

# Create a CSV table containing the number of mapped and paired reads
# for each sample
# Specialise this rule!!! using input functions
rule mappedAndPaired2csv:
    input:
        samples = 'samples/samples.tsv',
        HS678 = expand('output/HS678/alignment/target/celegans/ref/{{database}}/release{{databaseRelease}}/query/n2/{stage}/{bacterium}/repl{replicate}/hisat2/alnSortedByCoordStatsReadsMappedAndPaired.txt', stage = ['adults', 'embryos'], bacterium = ['BIGb446', 'BIGb468', 'hb101'], replicate = [1, 2, 3]),
        HS755N2 = expand('output/HS755/alignment/target/celegans/ref/{{database}}/release{{databaseRelease}}/query/n2/embryos/{bacterium}/repl{replicate}/hisat2/alnSortedByCoordStatsReadsMappedAndPaired.txt', bacterium = ['hb101', 'pa14', 'plum'], replicate = [1, 2, 3]),
        HS755NOB27 = expand('output/HS755/alignment/target/celegans/ref/{{database}}/release{{databaseRelease}}/query/nob27/embryos/hb101/repl{replicate}/hisat2/alnSortedByCoordStatsReadsMappedAndPaired.txt', replicate = [1, 2, 3]),
    output:
        csvTable = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/readsMappedPaired.csv',
    conda:
        'envs/conda/r-dplyr.yaml'
    script:
        'scripts/mappedAndPaired2csv.R'
# local
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS678/alignment/target/celegans/ref/ensembl/release96/query/readsMappedPaired.csv output/HS755/alignment/target/celegans/ref/ensembl/release96/query/readsMappedPaired.csv -n

rule quantify:
    input:
        index = 'output/transcriptome/celegans/ref/{database}/release{databaseRelease}/index/kallisto/transcriptome.idx',
        trimmed1 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/trimmed/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R1/L1L2.fq.gz',
        trimmed2 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/trimmed/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R2/L1L2.fq.gz',
    output:
        abundance = 'output/HS{HSrunNumber}/counts/{database}/release{databaseRelease}/kallisto/{wormStrain}/{stage}/{bacterium}/repl{replicate}/abundance.tsv',
    wildcard_constraints:
        database = '(ensembl)|(wormbase)',
        databaseRelease = '\d+',
    benchmark:
        'output/benchmark/quantify/HS{HSrunNumber}/{database}/release{databaseRelease}/{wormStrain}/{stage}{bacterium}/repl{replicate}.tsv'
    log:
        'output/HS{HSrunNumber}/counts/{database}/release{databaseRelease}/kallisto/{wormStrain}/{stage}/{bacterium}/repl{replicate}/quantify.log',
    conda:
        'envs/conda/kallisto0_45_1.yaml'
    threads: 8
    shell:
        'kallisto quant -i {input.index} -o output/HS{wildcards.HSrunNumber}/counts/{wildcards.database}/release{wildcards.databaseRelease}/kallisto/{wildcards.wormStrain}/{wildcards.stage}/{wildcards.bacterium}/repl{wildcards.replicate} --threads {threads} -b 100 {input.trimmed1} {input.trimmed2} > {log} 2>&1'
# local
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS755/counts/ensembl/release96/kallisto/n2/embryos/hb101/repl{1,2,3}/abundance.tsv -n
# SLURM Ensembl + Wormbase and HS678 + HS755
'''
snakemake --jobs 200 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason \
output/HS678/counts/ensembl/release{95,96}/kallisto/n2/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/abundance.tsv \
output/HS678/counts/wormbase/release{269,270}/kallisto/n2/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/abundance.tsv \
output/HS755/counts/ensembl/release{95,96}/kallisto/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/abundance.tsv \
output/HS755/counts/ensembl/release{95,96}/kallisto/nob27/embryos/hb101/repl{1,2,3}/abundance.tsv \
output/HS755/counts/wormbase/release{269,270}/kallisto/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/abundance.tsv \
output/HS755/counts/wormbase/release{269,270}/kallisto/nob27/embryos/hb101/repl{1,2,3}/abundance.tsv \
--cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} \
--job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} \
--ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} \
--mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n
'''

# Sort out this function and the next rule
# def findTPMs(wildcards):
#     HSrunNumber = int(wildcards[0])
#     if HSrunNumber == 678:
#     expand('output/HS{{HSrunNumber}}/counts/{{database}}/release{{databaseRelease}}/kallisto/n2/{stage}/{bacterium}/repl{replicate}/abundance.tsv', stage = ['adults', 'embryos'], bacterium = ['BIGb446', 'BIGb468', 'hb101'], replicate = [1, 2, 3]),
#     return

rule mergeTPMs:
    input:
        samples = 'samples/samples.tsv',
        tpm678 = expand('output/HS678/counts/{{database}}/release{{databaseRelease}}/kallisto/n2/{stage}/{bacterium}/repl{replicate}/abundance.tsv', stage = ['adults', 'embryos'], bacterium = ['BIGb446', 'BIGb468', 'hb101'], replicate = [1, 2, 3]),
        tpm755N2 = expand('output/HS755/counts/{{database}}/release{{databaseRelease}}/kallisto/n2/embryos/{bacterium}/repl{replicate}/abundance.tsv', stage = ['adults', 'embryos'], bacterium = ['hb101', 'pa14', 'plum'], replicate = [1, 2, 3]),
        tpm755NOB27 = expand('output/HS755/counts/{{database}}/release{{databaseRelease}}/kallisto/nob27/embryos/hb101/repl{replicate}/abundance.tsv', stage = ['adults', 'embryos'], bacterium = ['hb101', 'pa14', 'plum'], replicate = [1, 2, 3]),
    output:
        tpms = 'output/HS{HSrunNumber}/counts/{database}/release{databaseRelease}/kallisto/tpms.csv',
    conda:
        'envs/conda/r-dplyr.yaml'
    script:
        'scripts/mergeTPMsHS678.R'
# local
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS678/counts/ensembl/release96/kallisto/tpms.csv output/HS755/counts/ensembl/release96/kallisto/tpms.csv -n

# rule mergeTPMsHS755:
#     input:
#         n2tpm = expand('output/HS755/counts/{{database}}/release{{databaseRelease}}/kallisto/n2/embryos/{bacterium}/repl{replicate}/abundance.tsv', bacterium = ['hb101', 'pa14', 'plum'], replicate = [1, 2, 3]),
#         nob27tpm = expand('output/HS755/counts/{{database}}/release{{databaseRelease}}/kallisto/nob27/embryos/hb101/repl{replicate}/abundance.tsv', replicate = [1, 2, 3]),
#     output:
#         tpms = 'output/HS755/counts/{database}/release{databaseRelease}/kallisto/tpms.csv',
#     conda:
#         'envs/conda/r.yaml'
#     script:
#         'scripts/mergeTPMsHS755.R'
# local
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS755/counts/ensembl/release96/kallisto/tpms.csv -n

rule sortByName:
    input:
        bam = 'output/alignment/ensembl/release{ensemblRelease}/{stage}/{bact}/{repl}/Aligned.sortedByCoord.out.bam'
    output:
        bamSortedByName = 'output/alignment/ensembl/release{ensemblRelease}/{stage}/{bact}/{repl}/Aligned.sortedByName.out.bam'
    benchmark:
        'output/benchmark/sortByName/release{ensemblRelease}{stage}{bact}{repl}.tsv'
    log:
        'output/log/snakemake/sortByName/release{ensemblRelease}{stage}{bact}{repl}.txt'
    conda:
        'envs/conda/samtools1_9.yaml'
    threads: 6
    shell:
        'samtools sort -n --threads {threads} -o {output.bamSortedByName} {input.bam}'
# snakemake --jobs 200 output/alignment/ensembl/release95/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/Aligned.sortedByName.out.bam.bai --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}"

# indexing a BAMsortedByName does not work!
#rule indexBamSortedByName:
#    input:
#        bam = 'output/alignment/ensembl/release{ensemblRelease}/{stage}/{bact}/{repl}/Aligned.sortedByName.out.bam'
#    output:
#        index = 'output/alignment/ensembl/release{ensemblRelease}/{stage}/{bact}/{repl}/Aligned.sortedByName.out.bam.bai'
#    shell:
#        'samtools index {input.bam}'

rule featureCountsHisat2wormbase:
    input:
        genome = 'output/genome/celegans/ref/{database}/release{databaseRelease}/seq/c_elegans.PRJNA13758.genomic.fa',
        gtf = 'output/genome/celegans/ref/{database}/release{databaseRelease}/annotation/gff3/c_elegans.PRJNA13758.annotations.gff3',
        bamSortedByCoordN2 = expand('output/HS{{HSrunNumber}}/alignment/target/celegans/ref/{{database}}/release{{databaseRelease}}/query/n2/embryos/{bact}/repl{repl}/hisat2/alnSortedByCoord.bam', bact = ['hb101', 'pa14', 'plum'], repl = [1, 2, 3]),
        bamSortedByCoordNOB27 = expand('output/HS{{HSrunNumber}}/alignment/target/celegans/ref/{{database}}/release{{databaseRelease}}/query/nob27/embryos/hb101/repl{repl}/hisat2/alnSortedByCoord.bam',repl = [1, 2, 3]),
        index = expand('output/HS{{HSrunNumber}}/alignment/target/celegans/ref/{{database}}/release{{databaseRelease}}/query/n2/embryos/{bact}/repl{repl}/hisat2/alnSortedByCoord.bam.bai', bact = ['hb101', 'pa14', 'plum'], repl = [1, 2, 3]),
    output:
        counts = expand('output/HS{{HSrunNumber}}/counts/featurecounts/hisat2/celegans/{{database}}/release{{databaseRelease}}/minoverlap1FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt{extension}', extension = ['', '.summary', '.jcounts'])
    wildcard_constraints:
        database = '(ensembl)|(wormbase)',
        databaseRelease = '\d+',
    benchmark:
        'output/benchmark/featureCountsHisat2wormbase/HS{HSrunNumber}/{database}/release{databaseRelease}.tsv'
    log:
        'output/HS{HSrunNumber}/counts/featurecounts/hisat2/celegans/{database}/release{databaseRelease}/minoverlap1FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts_featureCountsHisat2all.log'
    conda:
        'envs/conda/subread1_6_3.yaml'
    threads: 8
    shell:
        'featureCounts -F GTF -t exon -g gene_id --minOverlap 1 -s 1 -J -G {input.genome} -p -B -P -d 50 -D 600 -T {threads} -a {input.gtf} -o {output.counts[0]} {input.bamSortedByCoordN2} > {log} 2>&1'
# -C makes it fail?
# {input.bamSortedByCoordNOB27}
# snakemake --cores 8 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS755/counts/featurecounts/hisat2/celegans/wormbase/release270/minoverlap1FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt -n
# SLURM
# snakemake --jobs 2000 --use-conda --resources downloads=1 --printshellcmds --keep-going --reason output/HS755/counts/featurecounts/hisat2/celegans/wormbase/release270/minoverlap1FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt{,.summary,.jcounts} output/HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt{,.summary,.jcounts} --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n

# featurecounts for HS678
rule featureCountsHisat2ensemblHS678:
    input:
        genome = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/seq/genome.fa',
        gtf = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/annotation/gtf/Caenorhabditis_elegans.WBcel235.gtf',
        bamSortedByCoord = expand('output/HS678/alignment/target/celegans/ref/ensembl/release{{ensemblRelease}}/query/n2/{stage}/{bact}/repl{repl}/hisat2/alnSortedByCoord.bam', stage = ['adults', 'embryos'], bact = ['hb101', 'BIGb446', 'BIGb468'], repl = [1, 2, 3]),
        index = expand('output/HS678/alignment/target/celegans/ref/ensembl/release{{ensemblRelease}}/query/n2/{stage}/{bact}/repl{repl}/hisat2/alnSortedByCoord.bam.bai', stage = ['adults', 'embryos'], bact = ['hb101', 'BIGb446', 'BIGb468'], repl = [1, 2, 3]),
    output:
        counts = expand('output/HS678/counts/featurecounts/hisat2/celegans/ensembl/release{{ensemblRelease}}/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt{extension}', extension = ['', '.summary', '.jcounts'])
    wildcard_constraints:
        ensemblRelease = '\d*'
    benchmark:
        'output/benchmark/featureCountsHisat2ensemblStranded/release{ensemblRelease}.tsv'
    log:
        'output/HS678/counts/featurecounts/hisat2/celegans/ensembl/release{ensemblRelease}/featureCountsHisat2ensemblStranded.log'
    conda:
        'envs/conda/subread1_6_3.yaml'
    threads: 8
    shell:
        'featureCounts -F GTF -t exon -g gene_id --minOverlap 1 -J -G {input.genome} -p -B -P -d 50 -D 600 -T {threads} -a {input.gtf} -o {output.counts[0]} {input.bamSortedByCoord} > {log} 2>&1'
# -C
# SLURM
# snakemake --jobs 200 --use-conda --resources downloads=1 --printshellcmds --keep-going --reason output/HS678/counts/featurecounts/hisat2/celegans/ensembl/release{95,96}/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt{,.summary,.jcounts} --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n

rule featureCountsHisat2ensemblHS755:
    input:
        genome = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/seq/genome.fa',
        gtf = 'output/genome/celegans/ref/ensembl/release{ensemblRelease}/annotation/gtf/Caenorhabditis_elegans.WBcel235.gtf',
        bamSortedByCoordN2 = expand('output/HS755/alignment/target/celegans/ref/ensembl/release{{ensemblRelease}}/query/n2/embryos/{bact}/repl{repl}/hisat2/alnSortedByCoord.bam', bact = ['hb101', 'pa14', 'plum'], repl = [1, 2, 3]),
        bamSortedByCoordNOB27 = expand('output/HS755/alignment/target/celegans/ref/ensembl/release{{ensemblRelease}}/query/nob27/embryos/hb101/repl{repl}/hisat2/alnSortedByCoord.bam',repl = [1, 2, 3]),
        index = expand('output/HS755/alignment/target/celegans/ref/ensembl/release{{ensemblRelease}}/query/n2/embryos/{bact}/repl{repl}/hisat2/alnSortedByCoord.bam.bai', bact = ['hb101', 'pa14', 'plum'], repl = [1, 2, 3]),
    output:
        counts = expand('output/HS755/counts/featurecounts/hisat2/celegans/ensembl/release{{ensemblRelease}}/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt{extension}', extension = ['', '.summary', '.jcounts'])
    wildcard_constraints:
        ensemblRelease = '\d*'
    benchmark:
        'output/benchmark/featureCountsHisat2ensembl/release{ensemblRelease}.tsv'
    log:
        'output/HS755/counts/featurecounts/hisat2/celegans/ensembl/release{ensemblRelease}/minoverlap1FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts_featureCountsHisat2all.log'
    conda:
        'envs/conda/subread1_6_3.yaml'
    threads: 8
    shell:
        'featureCounts -F GTF -t exon -g gene_id --minOverlap 1 -s 0 -J -G {input.genome} -p -B -P -d 50 -D 600 -T {threads} -a {input.gtf} -o {output.counts[0]} {input.bamSortedByCoordN2} {input.bamSortedByCoordNOB27} > {log} 2>&1'
# -C makes it fail?
# {input.bamSortedByCoordNOB27}
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS755/counts/featurecounts/hisat2/celegans/ensembl/release270/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt -n
# SLURM
# snakemake --jobs 2000 --use-conda --resources downloads=1 --printshellcmds --keep-going --reason output/HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt{,.summary,.jcounts} --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n

rule featureCountsHisat2ensemblHS678HS755:
    input:
        genome = 'output/genome/celegans/ref/{database}/release{ensemblRelease}/seq/genome.fa',
        gtf = 'output/genome/celegans/ref/{database}/release{ensemblRelease}/annotation/gtf/Caenorhabditis_elegans.WBcel235.gtf',
        bamSortedByCoordHS678 = expand('output/HS678/alignment/target/celegans/ref/{{database}}/release{{ensemblRelease}}/query/n2/{stage}/{bact}/repl{repl}/hisat2/alnSortedByCoord.bam', stage = ['adults', 'embryos'], bact = ['hb101', 'BIGb446', 'BIGb468'], repl = [1, 2, 3]),
        indexHS678 = expand('output/HS678/alignment/target/celegans/ref/{{database}}/release{{ensemblRelease}}/query/n2/{stage}/{bact}/repl{repl}/hisat2/alnSortedByCoord.bam.bai', stage = ['adults', 'embryos'], bact = ['hb101', 'BIGb446', 'BIGb468'], repl = [1, 2, 3]),
        bamSortedByCoordN2 = expand('output/HS755/alignment/target/celegans/ref/{{database}}/release{{ensemblRelease}}/query/n2/embryos/{bact}/repl{repl}/hisat2/alnSortedByCoord.bam', bact = ['hb101', 'pa14', 'plum'], repl = [1, 2, 3]),
        bamSortedByCoordNOB27 = expand('output/HS755/alignment/target/celegans/ref/{{database}}/release{{ensemblRelease}}/query/nob27/embryos/hb101/repl{repl}/hisat2/alnSortedByCoord.bam',repl = [1, 2, 3]),
        indexHS755 = expand('output/HS755/alignment/target/celegans/ref/{{database}}/release{{ensemblRelease}}/query/n2/embryos/{bact}/repl{repl}/hisat2/alnSortedByCoord.bam.bai', bact = ['hb101', 'pa14', 'plum'], repl = [1, 2, 3]),
    output:
        counts = expand('output/HS678HS755/counts/featurecounts/hisat2/celegans/{{database}}/release{{ensemblRelease}}/minoverlap1{{strandedness}}FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt{extension}', extension = ['', '.summary', '.jcounts'])
    wildcard_constraints:
        database = 'ensembl',
        ensemblRelease = '\d*',
    benchmark:
        'output/benchmark/featureCountsHisat2ensemblHS678HS755/{database}/release{ensemblRelease}/{strandedness}.tsv'
    log:
        'output/HS678HS755/counts/featurecounts/hisat2/celegans/{database}/release{ensemblRelease}/minoverlap1{strandedness}FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts_featureCountsHisat2all.log'
    conda:
        'envs/conda/subread1_6_3.yaml'
    threads: 8
    shell:
        'featureCounts -F GTF -t exon -g gene_id --minOverlap 1 -s 0 -J -G {input.genome} -p -B -P -d 50 -D 600 -T {threads} -a {input.gtf} -o {output.counts[0]} {input.bamSortedByCoordHS678} {input.bamSortedByCoordN2} {input.bamSortedByCoordNOB27} > {log} 2>&1'
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS678HS755/counts/featurecounts/hisat2/celegans/ensembl/release270/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt -n
# SLURM
# snakemake --jobs 2000 --use-conda --resources downloads=1 --printshellcmds --keep-going --reason output/HS678HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt{,.summary,.jcounts} --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n

# Rename the names of the columns of the counts table to the sample name
# (variable1_variable2_replicate)
rule renameFeaturecounts:
    input:
        conditions = 'samples/conditions.tsv',
        counts = 'output/HS{HSrunNumber}/counts/featurecounts/hisat2/celegans/{database}/release{databaseRelease}/minoverlap1{strandedness}FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCounts.txt',
    output:
        counts = 'output/HS{HSrunNumber}/counts/featurecounts/hisat2/celegans/{database}/release{databaseRelease}/minoverlap1{strandedness}FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsRenamed.txt',
    benchmark:
        'output/benchmark/renameFeaturecounts/HS{HSrunNumber}/{database}/{databaseRelease}/{strandedness}.tsv'
    log:
        'output/HS{HSrunNumber}/counts/featurecounts/hisat2/celegans/{database}/release{databaseRelease}/{strandedness}_renameFeaturecounts.log',
    conda:
        'envs/conda/r-dplyr_r-stringr.yaml'
    script:
        'scripts/renameFeaturecounts.R'
# local
# snakemake --cores 7 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS678/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsRenamed.txt output/HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsRenamed.txt -n
# SLURM
# snakemake --jobs 2000 --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS678/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsRenamed.txt output/HS755/counts/featurecounts/hisat2/celegans/ensembl/release96/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsRenamed.txt --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n


# Calculate FPM for each gene in each sample
# FPM = Fragment Per Million = number of mapped paired-end fragment hitting an exon in a gene divided by the total
# number of mapped paired-end fragments
rule calculateFPKMs:
    input:
        readsMappedPaired = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/{database}/release{databaseRelease}/query/readsMappedPaired.csv',
        counts = 'output/HS{HSrunNumber}/counts/featurecounts/hisat2/celegans/{database}/release{databaseRelease}/minoverlap1{strandedness}FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersCountsRenamed.txt',
    output:
        fpkm = 'output/HS{HSrunNumber}/counts/featurecounts/hisat2/celegans/{database}/release{databaseRelease}/minoverlap1{strandedness}FragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersFPKMs.txt',
    wildcard_constraints:
        strandedness = '(stranded)|(unstranded)'
    benchmark:
        'output/benchmark/calculateFPKMs/HS{HSrunNumber}/{database}/{databaseRelease}/{strandedness}.tsv'
    log:
        'output/HS{HSrunNumber}/counts/featurecounts/hisat2/celegans/{database}/release{databaseRelease}/{strandedness}_calculateFPKMs.log',
    conda:
        'envs/conda/r-dplyr_r-stringr.yaml'
    script:
        'scripts/calculateFPKMs.R'
# local
# snakemake --cores 8  --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS{678,755}/counts/featurecounts/hisat2/celegans/ensembl/release{95,96}/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersFPKMs.txt -n
# SLURM
# snakemake --jobs 2000  --resources downloads=1 --use-conda --printshellcmds --keep-going --reason output/HS{678,755}/counts/featurecounts/hisat2/celegans/ensembl/release{95,96}/minoverlap1unstrandedFragmentsBothReadsAlignedMinFragLength50MaxFragLength600noMultimappersFPKMs.txt --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}" -n

rule featureCountsAll:
    input:
        genome = 'output/genome/ensembl/release{ensemblRelease}/seq/genome.fa.gz',
        gtf = 'output/genome/ensembl/release{ensemblRelease}/annotation/Caenorhabditis_elegans.WBcel235.gtf',
        bamSortedByCoord = expand('output/alignment/ensembl/release{{ensemblRelease}}/{stage}/{bact}/{repl}/Aligned.sortedByCoord.out.bam', stage = stages, bact = bacteria, repl = replicates),
        index = expand('output/alignment/ensembl/release{{ensemblRelease}}/{stage}/{bact}/{repl}/Aligned.sortedByCoord.out.bam.bai', stage = stages, bact = bacteria, repl = replicates),
    output:
        counts = report('output/counts/featurecounts/ensembl/release{ensemblRelease}/counts.txt')
    wildcard_constraints:
        ensemblRelease = '\d*'
    benchmark:
        'output/benchmark/featureCountsAll/release{ensemblRelease}.tsv'
    log:
        'output/log/featureCountsAll/release{ensemblRelease}.txt'
    conda:
        'envs/conda/subread.yaml'
    threads: coresPerNode
    shell:
        'featureCounts -F GTF -t exon -g gene_id --minOverlap 1 -s 1 -J -G {input.genome} -p -T {threads} -a {input.gtf} -o {output.counts} {input.bamSortedByCoord}'

# -M Multi-mapping reads will also be counted. For a multi-mapping read, all its reported
# alignments will be counted. The 'NH' tag in BAM/SAM input is used to detect multi-mapping reads.
# --fraction Assign fractional counts to features. This option must
#  be used together with '-M' or '-O' or both. When '-M' is
#  specified, each reported alignment from a multi-mapping
#  read (identified via 'NH' tag) will carry a fractional
#  count of 1/x, instead of 1 (one), where x is the total
#  number of alignments reported for the same read. When '-O'
#  is specified, each overlapping feature will receive a
#  fractional count of 1/y, where y is the total number of
#  features overlapping with the read. When both '-M' and
#  '-O' are specified, each alignment will carry a fractional
#  count of 1/(x*y).
rule featureCountsAllMultimappersMultioverlappers:
    input:
        genome = 'output/genome/ensembl/release{ensemblRelease}/seq/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa',
        gtf = 'output/genome/ensembl/release{ensemblRelease}/annotation/Caenorhabditis_elegans.WBcel235.gtf',
        bamSortedByCoord = expand('output/alignment/ensembl/release{{ensemblRelease}}/{stage}/{bact}/{repl}/Aligned.sortedByCoord.out.bam', stage = stages, bact = bacteria, repl = replicates),
        index = expand('output/alignment/ensembl/release{{ensemblRelease}}/{stage}/{bact}/{repl}/Aligned.sortedByCoord.out.bam.bai', stage = stages, bact = bacteria, repl = replicates),
    output:
        counts = report('output/counts/featurecounts/ensembl/release{ensemblRelease}/countsMultimappersMultioverlappers.txt')
    wildcard_constraints:
        ensemblRelease = '\d*'
    benchmark:
        'output/benchmark/featureCountsAll/release{ensemblRelease}.tsv'
    log:
        'output/log/featureCountsAll/release{ensemblRelease}.txt'
    conda:
        'envs/conda/subread.yaml'
    threads: coresPerNode
    shell:
        'featureCounts -F GTF -t exon -g gene_id --minOverlap 1 -M --fraction -s 1 -J -G {input.genome} -p -T {threads} -a {input.gtf} -o {output.counts} {input.bamSortedByCoord}'

rule htseqCount:
    input:
        gtf = 'output/genome/ensembl/release{ensemblRelease}/annotation/Caenorhabditis_elegans.WBcel235.gtf',
        bam = 'output/alignment/ensembl/release{ensemblRelease}/{stage}/{bact}/{repl}/Aligned.sortedByName.out.bam',
    output:
        counts = 'output/counts/ensembl/release{ensemblRelease}/htseq/{stage}/{bact}/{repl}/counts.txt'
    benchmark:
        'output/benchmark/htseqCount/release{ensemblRelease}{stage}{bact}{repl}.tsv'
    log:
        'output/log/snakemake/htseqCount/release{ensemblRelease}{stage}{bact}{repl}.txt'
    conda:
        'envs/conda/htseq.yaml'
    shell:
        'htseq-count --mode=intersection-nonempty --stranded=yes --type=exon --idattr=gene_id --order=name --format=bam {input.bam} {input.gtf} > {output.counts}'

rule embryosDds2:
    input:
        counts = 'output/counts/featurecounts/ensembl/release{ensemblRelease}/counts.txt'
    output:
        dds2 = 'output/robjects/featurecounts/ensembl/release{ensemblRelease}/dds2.Rda',
        #heatmap = 'output/images/heatmap/ensembl/release{ensemblRelease}/adultCorrelations.pdf',
        #tree = 'output/images/tree/ensembl/release{ensemblRelease}/embryosHierarchicalClustering.pdf',
        #scatter = 'output/images/scatter/ensembl/release{ensemblRelease}/embryosPCArnaExprStarDESeq2.pdf',
        #dge = 'output/DGE/ensembl/release{ensemblRelease}/embryosPC1.csv'
    benchmark:
        'output/benchmark/DGEembryos/ensembl/release{ensemblRelease}.tsv'
    log:
        'output/log/snakemake/DGEembryos/ensembl/release{ensemblRelease}.txt'
    conda:
        'envs/conda/bioconductor-deseq2_r-stringr.yaml'
    script:
        'scripts/04dgeEmbryos.R'

rule irods2lustre:
    output:
        'output/reads/dna/pacbio/all/a.subreads.bam',
    shell:
        'bash scripts/01download.sh'
# snakemake output/reads/dna/pacbio/all/a.subreads.bam

rule bam2fqPacbio:
    input:
        bam = 'output/reads/dna/pacbio/{pathogenORall}/a.subreads.bam'
    output:
        fq = 'output/reads/dna/pacbio/{pathogenORall}/a.subreads.fq'
    log:
        'output/log/conda/bam2fqPacbio{pathogenORall}.txt'
    conda:
        'envs/conda/bamtools.yaml'
    shell:
        'bamtools convert -format fastq -in {input.bam} -out {output.fq}'
# snakemake --jobs 4 --use-conda output/reads/dna/pacbio/BIGb248/a.subreads.fq --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}"

rule fq2fa:
    input:
        fq = 'output/reads/dna/pacbio/{pathogenORall}/a.subreads.fq'
    output:
        fa = 'output/reads/dna/pacbio/{pathogenORall}/a.subreads.fa'
    log:
        'output/log/conda/fq2fasta{pathogenORall}.txt'
    conda:
        'envs/conda/seqtk.yaml'
    shell:
        'seqtk seq -a {input.fq} > {output.fa}'
# snakemake --use-conda output/reads/dna/pacbio/all/a.subreads.fa

# conda install -c bioconda lima

rule assemble:
    input:
        reads = 'output/reads/dna/pacbio/{pathogen}/a.subreads.fq'
    output:
        assembly = 'output/assembly/canu/{pathogen}/assembly.contigs.fasta'
    log:
        'output/log/conda/assemble{pathogen}.txt'
    threads: 42 # recommended by LSF
    conda:
        'envs/conda/canu.yaml'
    shell:
        'canu -p {wildcards.pathogen} -d output/assembly/canu/{wildcards.pathogen} genomeSize=6m useGrid=false -pacbio-raw {input.reads}'

# rule annotateGenomes:
#     input:
#         genome = 'input/assembly/fasta/BIGb248.fasta'
#     output:
#     log:
#     threads: 8
#     conda:
#         'envs/conda/prokka.yaml'
#     shell:
#         'prokka {input.genome} --outdir output/assembly/hgap4/annotation/prokka/BIGb248 --rfam --force'

# on laptop:
# canu -p BIGb248 -d output/assembly/canu/BIGb248 genomeSize=6m -pacbio-raw output/reads/dna/pacbio/BIGb248/a.subreads.fq
# snakemake --cores 7 --use-conda --printshellcmds output/assembly/canu/{BIGb248,BIGb306,BIGb446,BIGb468}/assembly.contigs.fasta
# on cluster:
# bsub -q basement -E test -e /nfs/users/nfs_c/cr7 -n20 -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M100000 -o canu.o -e canu.e -J canu 'bash canu.sh'
# snakemake --jobs 4 --use-conda --printshellcmds output/assembly/canu/BIGb248/assembly.contigs.fasta --cluster-config slurm.json --cluster "sbatch -A {cluster.A} -p {cluster.p} -t {cluster.time} --job-name {cluster.jobName} --output {cluster.output} --error {cluster.error} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mail-type {cluster.mailType} --mail-user {cluster.mailUser}"

# conda create --name gff3toembl --channel bioconda --yes genometools gff3toembl
# conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/gff3toembl.yaml
# gff3_to_embl celegans 6239 5569 pathogens output/assembly/hgap4/annotation/prokka/BIGb248/PROKKA_02072019.gff
