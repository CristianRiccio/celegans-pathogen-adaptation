rule hisat2wormbase:
    input:
        index = expand('output/genome/celegans/ref/wormbase/release{{wormbaseRelease}}/index/hisat2/index.{number}.ht2', number = ['1', '2', '3', '4', '5', '6', '7', '8']),
        reads1 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/trimmed/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R1/L1L2.fq.gz',
        reads2 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/trimmed/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R2/L1L2.fq.gz',
    output:
        sam = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/wormbase/release{wormbaseRelease}/query/{wormStrain}/{stage}/{bacterium}/repl{replicate}/hisat2/aln.sam'
    wildcard_constraints:
        wormbaseRelease = '\d+'
    benchmark:
        'output/benchmark/hisat2wormbase/HS{HSrunNumber}/release{wormbaseRelease}/{wormStrain}/{stage}/{bacterium}/{replicate}.tsv'
    log:
        'output/HS{HSrunNumber}/alignment/target/celegans/ref/ensembl/release{wormbaseRelease}/query/{wormStrain}/{stage}/{bacterium}/{replicate}/hisat2/hisat2wormbase.log'
    threads: 8
    conda:
        '../envs/conda/hisat2v2_1_0.yaml'
    shell:
        'hisat2 -k 5000 --threads {threads} -1 {input.reads1} -2 {input.reads2} -x output/genome/celegans/ref/wormbase/release{wildcards.wormbaseRelease}/index/hisat2/index -S {output.sam} > {log} 2>&1'
# local
# snakemake --cores 8 output/HS755/alignment/target/celegans/ref/wormbase/release270/query/nob27/embryos/hb101/repl{1,2,3}/hisat2/aln.sam output/HS755/alignment/target/celegans/ref/wormbase/release270/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/aln.sam -n
'''
snakemake --profile slurm output/HS678/alignment/target/celegans/ref/wormbase/release270/query/n2/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/aln.sam \
output/HS755/alignment/target/celegans/ref/wormbase/release270/query/nob27/embryos/hb101/repl{1,2,3}/hisat2/aln.sam \
output/HS755/alignment/target/celegans/ref/wormbase/release270/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/aln.sam \
-n
'''

# snakemake --cores 8 output/alignment/target/celegans/ref/ensembl/release95/query/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/aln.sam
# snakemake --profile slurm output/alignment/target/celegans/ref/ensembl/release95/query/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/hisat2/aln.sam -n

rule hisat2ensembl:
    input:
        index = expand('output/genome/celegans/ref/ensembl/release{{ensemblRelease}}/index/hisat2/index.{number}.ht2', number = ['1', '2', '3', '4', '5', '6', '7', '8']),
        reads1 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/trimmed/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R1/L1L2.fq.gz',
        reads2 = 'output/HS{HSrunNumber}/reads/rna/illumina/celegans/trimmed/{wormStrain}/{stage}/{bacterium}/repl{replicate}/R2/L1L2.fq.gz',
    output:
        sam = 'output/HS{HSrunNumber}/alignment/target/celegans/ref/ensembl/release{ensemblRelease}/query/{wormStrain}/{stage}/{bacterium}/repl{replicate}/hisat2/aln.sam'
    wildcard_constraints:
        ensemblRelease = '\d+'
    benchmark:
        'output/benchmark/hisat2ensembl/HS{HSrunNumber}/release{ensemblRelease}/{wormStrain}/{stage}/{bacterium}/repl{replicate}.tsv'
    log:
        'output/HS{HSrunNumber}/alignment/target/celegans/ref/ensembl/release{ensemblRelease}/query/{wormStrain}/{stage}/{bacterium}/repl{replicate}/hisat2/hisat2ensembl.log'
    threads: 8
    conda:
        '../envs/conda/hisat2v2_1_0.yaml'
    shell:
        'hisat2 -k 5000 --threads {threads} -1 {input.reads1} -2 {input.reads2} -x output/genome/celegans/ref/ensembl/release{wildcards.ensemblRelease}/index/hisat2/index -S {output.sam} > {log} 2>&1'
# snakemake --cores 8 output/HS755/alignment/target/celegans/ref/ensembl/release96/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/aln.sam
'''
snakemake --profile slurm output/HS678/alignment/target/celegans/ref/ensembl/release{95,96}/query/n2/{adults,embryos}/{hb101,BIGb446,BIGb468}/repl{1,2,3}/hisat2/aln.sam \
output/HS755/alignment/target/celegans/ref/ensembl/release96/query/n2/embryos/{hb101,pa14,plum}/repl{1,2,3}/hisat2/aln.sam \
output/HS755/alignment/target/celegans/ref/ensembl/release96/query/nob27/embryos/hb101/repl{1,2,3}/hisat2/aln.sam \
-n
'''


rule STARensembl:
    input:
        index = expand('output/genome/ensembl/release{{ensemblRelease}}/index/{extension}', extension = ['Genome', 'SAindex', 'SA', 'chrName.txt', 'chrStart.txt', 'chrLength.txt', 'chrNameLength.txt', 'genomeParameters.txt']),
        gtf = 'output/genome/ensembl/release{ensemblRelease}/annotation/Caenorhabditis_elegans.WBcel235.gtf',
        reads1 = 'output/reads/trimmed/{stage}/{bact}/{repl}/R1.fq.gz',
        reads2 = 'output/reads/trimmed/{stage}/{bact}/{repl}/R2.fq.gz'
    output:
        bamSortedByCoord = 'output/alignment/ensembl/release{ensemblRelease}/{stage}/{bact}/{repl}/Aligned.sortedByCoord.out.bam'
    benchmark:
        'output/benchmark/STAR/release{ensemblRelease}{stage}{bact}{repl}.tsv'
    log:
        'output/log/snakemake/STAR/release{ensemblRelease}{stage}{bact}{repl}.txt'
    threads: 8
    conda:
        '../envs/conda/star2_7_0.yaml'
    shell:
        'STAR --runMode alignReads --runThreadN {threads}  --outSAMtype BAM  SortedByCoordinate Unsorted --sjdbOverhang 100  --genomeDir output/genome/ensembl/release{wildcards.ensemblRelease}/index --readFilesIn {input.reads1} {input.reads2} --readFilesCommand gunzip -c --outFileNamePrefix output/alignment/ensembl/release{wildcards.ensemblRelease}/{wildcards.stage}/{wildcards.bact}/{wildcards.repl}/ --sjdbGTFfile {input.gtf}'
# snakemake --cores 7 output/alignment/ensembl/release96/adults/hb101/repl1/Aligned.sortedByCoord.out.bam -n

rule STARwormbase:
    input:
        index = expand('output/genome/celegans/ref/wormbase/release{{wormbaseRelease}}/index/withGFF3/{extension}', extension = ['Genome', 'SAindex', 'SA', 'chrName.txt', 'chrStart.txt', 'chrLength.txt', 'chrNameLength.txt', 'genomeParameters.txt']),
        gff3 = 'output/genome/celegans/ref/wormbase/release{wormbaseRelease}/annotation/gff3/c_elegans.PRJNA13758.annotations.gff3',
        reads1 = 'output/reads/trimmed/{stage}/{bact}/{repl}/R1.fq.gz',
        reads2 = 'output/reads/trimmed/{stage}/{bact}/{repl}/R2.fq.gz'
    output:
        bamSortedByCoord = 'output/alignment/target/celegans/ref/wormbase/release{wormbaseRelease}/query/{stage}/{bact}/{repl}/Aligned.sortedByCoord.out.bam'
    benchmark:
        'output/benchmark/STARwormbase/release{wormbaseRelease}{stage}{bact}{repl}.tsv'
    log:
        'output/log/snakemake/STARwormbase/release{wormbaseRelease}{stage}{bact}{repl}.txt'
    threads: 8
    conda:
        '../envs/conda/star2_7_0.yaml'
    shell:
        'STAR --runMode alignReads --runThreadN {threads}  --outSAMtype BAM  SortedByCoordinate Unsorted --sjdbOverhang 49  --genomeDir output/genome/celegans/ref/wormbase/release{wildcards.wormbaseRelease}/index/withGFF3 --readFilesIn {input.reads1} {input.reads2} --readFilesCommand gunzip -c --outFileNamePrefix output/alignment/target/celegans/ref/wormbase//release{wildcards.wormbaseRelease}/query/{wildcards.stage}/{wildcards.bact}/{wildcards.repl}/ --sjdbGTFfile {input.gff3}'
# snakemake --profile slurm output/alignment/target/celegans/ref/wormbase/release269/query/{adults,embryos}/{BIGb446,BIGb468,hb101}/repl{1,2,3}/Aligned.sortedByCoord.out.bam -n
