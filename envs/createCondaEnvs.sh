#### Create conda environments

# nick
conda create --name nick --channel conda-forge --channel bioconda --yes git=2.21.0 r=3.5.1 snakemake=5.4.0
source activate nick
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/nick.yaml

source deactivate

# blat
conda create --name blat --channel bioconda --yes blat=36
source activate blat
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/blat.yaml

source deactivate

# bcftools
conda create --name bcftools --channel conda-forge --channel bioconda --yes bcftools=1.9
source activate bcftools
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bcftools.yaml

source deactivate

# bedtools
conda create --name bedtools --channel conda-forge --channel bioconda --yes bedtools=2.27.1
source activate bedtools
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bedtools.yaml

source deactivate

bioconductor-biomart_bioconductor-deseq2_r-stringr_r-tidyr
conda create --name bioconductor-biomart_bioconductor-deseq2_r-stringr_r-tidyr --channel conda-forge --yes bioconductor-biomart=2.38.0 bioconductor-deseq2=1.22.1 r-stringr=1.3.1 r-tidyr=0.8.3
source activate bioconductor-biomart_bioconductor-deseq2_r-stringr_r-tidyr
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-biomart_bioconductor-deseq2_r-stringr_r-tidyr.yaml

source deactivate

# bioconductor-deseq2
conda create --name bioconductor-deseq2 --channel conda-forge --channel bioconda --yes bioconductor-deseq2=1.22.1
source activate bioconductor-deseq2
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-deseq2.yaml

source deactivate

# bioconductor-deseq2_r-stringr
conda create --name bioconductor-deseq2_r-stringr --channel conda-forge --channel bioconda --channel r --yes bioconductor-deseq2=1.22.1 r=3.5.1 r-stringr=1.3.1
source activate bioconductor-deseq2_r-stringr
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-deseq2_r-stringr.yaml

source deactivate

# bioconductor-deseq2_r-dplyr_r-stringr
conda create --name bioconductor-deseq2_r-dplyr_r-stringr --channel conda-forge --channel bioconda --yes bioconductor-deseq2=1.22.1 r=3.5.1 r-dplyr=0.8.0 r-stringr=1.3.1
source activate bioconductor-deseq2_r-dplyr_r-stringr
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-deseq2_r-dplyr_r-stringr.yaml

source deactivate

# bioconductor-deseq2_r-plotly_r-stringr
conda create --name bioconductor-deseq2_r-plotly_r-stringr --channel conda-forge --channel bioconda --yes bioconductor-deseq2=1.22.1 r=3.5.1 r-plotly=4.8.0 r-stringr=1.3.1
source activate bioconductor-deseq2_r-plotly_r-stringr
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-deseq2_r-plotly_r-stringr.yaml

source deactivate

# bioconductor-deseq2_r-dplyr_r-pheatmap_r-stringr
conda create --name bioconductor-deseq2_r-dplyr_r-pheatmap_r-stringr --channel conda-forge --channel bioconda --yes bioconductor-deseq2=1.22.1 r=3.5.1 r-dplyr=0.8.0 r-pheatmap=1.0.12 r-stringr=1.3.1
source activate bioconductor-deseq2_r-dplyr_r-pheatmap_r-stringr
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-deseq2_r-dplyr_r-pheatmap_r-stringr.yaml

source deactivate

# bioconductor-rtracklayer
conda create --name bioconductor-rtracklayer --channel conda-forge --channel bioconda --channel r --yes bioconductor-rtracklayer=1.42.1 r=3.5.1
source activate bioconductor-rtracklayer
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-rtracklayer.yaml

source deactivate

# bioconductor-rsamtools_r-seqinr_r-tidyr
conda create --name bioconductor-rsamtools_r-seqinr_r-tidyr --channel conda-forge --channel bioconda --yes bioconductor-rsamtools=1.34.0 r=3.5.1 r-seqinr=3.4_5 r-tidyr=0.8.2
conda activate bioconductor-rsamtools_r-seqinr_r-tidyr
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-rsamtools_r-seqinr_r-tidyr.yaml

source deactivate

# bioconductor-rtracklayer_r-dplyr_r-ggplot2
conda create --name bioconductor-rtracklayer_r-dplyr_r-ggplot2 --channel conda-forge --channel bioconda --yes r=3.5.1 bioconductor-rtracklayer=1.42.1 r-dplyr=0.8.0 r-ggplot2=3.1.0
source activate bioconductor-rtracklayer_r-dplyr_r-ggplot2
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bioconductor-rtracklayer_r-dplyr_r-ggplot2.yaml

source deactivate

# bwa
conda create --name bwa --channel bioconda --yes bwa=0.7.17
source activate bwa
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bwa.yaml

source deactivate

# canu, version on Linux ahead of Mac
conda create --name canu --channel bioconda --yes canu=1.8
source activate canu
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/canu.yaml

source deactivate

# csvkit, openpyxl 2.6.0 breaks csvkit
conda create --name csvkit --channel conda-forge --yes csvkit=1.0.3 openpyxl=2.5.0
source activate csvkit
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/csvkit.yaml

source deactivate

# cufflinks2_2_1
conda create --name cufflinks2_2_1 --channel conda-forge --channel bioconda --channel r --yes cufflinks=2.2.1
source activate cufflinks2_2_1
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/cufflinks2_2_1.yaml

source deactivate

# cutadapt1_18
conda create --name cutadapt1_18 --channel bioconda --yes cutadapt=1.18
source activate cutadapt1_18
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/cutadapt1_18.yaml

source deactivate

# fastqc0_11_8
conda create --name fastqc0_11_8 --channel conda-forge --channel bioconda --yes fastqc=0.11.8
source activate fastqc0_11_8
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/fastqc0_11_8.yaml

source deactivate

# fastx_toolkit
conda create --name fastx_toolkit --channel bioconda --yes fastx_toolkit=0.0.14
source activate fastx_toolkit
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/fastx_toolkit.yaml

source deactivate

# hisat2v2_1_0
conda create --name hisat2v2_1_0 --channel bioconda --channel conda-forge --yes hisat2=2.1.0
source activate hisat2v2_1_0
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/hisat2v2_1_0.yaml

source deactivate

# htseq
conda create --name htseq --channel conda-forge --channel bioconda --yes htseq=0.11.2
conda activate htseq
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/htseq.yaml

source deactivate


# kallisto0_45_1
conda create --name kallisto0_45_1 --channel conda-forge --channel bioconda --yes kallisto=0.45.1
source activate kallisto0_45_1
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/kallisto0_45_1.yaml

source deactivate

# minimap2
conda create --name minimap2 --channel bioconda --yes minimap2=2.15
source activate minimap2
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/minimap2.yaml

source deactivate

# multiqc1_7
conda create --name multiqc1_7 --channel conda-forge --channel bioconda --yes multiqc=1.7
source activate multiqc1_7
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/multiqc1_7.yaml

source deactivate

# picard
conda create --name picard --channel bioconda --yes picard=2.18.26
source activate picard
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/picard.yaml

source deactivate

# python3_7_3
conda create --name python3_7_3 --channel conda-forge --yes python=3.7.3
source activate python3_7_3
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/python3_7_3.yaml

source deactivate

# plotQC_N2
conda create --name plotQC_N2 --channel conda-forge --yes r r-dplyr r-ggplot2
source activate plotQC_N2
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/plotQC_N2.yaml

source deactivate

# r3_5_1
conda create --name r --channel conda-forge --yes r=3.5.1
source activate r3_5_1
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r3_5_1.yaml

# r-all
conda create --name r-all --channel conda-forge --yes  bioconductor-deseq2=1.22.1 bioconductor-rtracklayer=1.42.1 bioconductor-rsamtools=1.34.0 r=3.5.1 r-dplyr=0.8.1 r-ggplot2=3.1.0 r-plotly=4.8.0 r-seqinr=3.4_5 r-stringr=1.4.0 r-tidyr=0.8.3 r-yaml=2.2.0
source activate r-all
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-all.yaml

source deactivate

# r-dplyr
conda create --name r-dplyr --channel conda-forge --yes r=3.5.1 r-dplyr=0.8.1
source activate r-dplyr
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-dplyr.yaml

source deactivate

# r-dplyr_r-stringr
conda create --name r-dplyr_r-stringr --channel conda-forge --yes r=3.5.1 r-dplyr=0.8.1  r-stringr=1.4.0
source activate r-dplyr_r-stringr
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-dplyr_r-stringr.yaml

source deactivate

# r-ggplot2
conda create --name r-ggplot2 --channel conda-forge --yes r=3.5.1 r-ggplot2=3.1.0
source activate r-ggplot2
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-ggplot2.yaml

source deactivate

# r-ggplot2_r-plotly_r-stringr_r-tidyr
conda create --name r-ggplot2_r-plotly_r-stringr_r-tidyr --channel conda-forge --yes r=3.5.1 r-ggplot2=3.1. r-plotly=4.8.0 r-stringr=1.3.1 r-tidyr=0.8.3
source activate r-ggplot2_r-plotly_r-stringr_r-tidyr
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-ggplot2_r-plotly_r-stringr_r-tidyr.yaml

source deactivate

# r-plotly_r-tidyverse
conda create --name r-ggplot2 --channel conda-forge --yes r=3.5.1 r-plotly=4.8.0 r-tidyverse=1.2.1
source activate r-plotly_r-tidyverse
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-plotly_r-tidyverse.yaml

source deactivate

# r-seqinr
conda create --name r-seqinr --channel conda-forge --yes r=3.5.1 r-seqinr=3.4_5
source activate r-seqinr
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-seqinr.yaml

source deactivate

# r-stringr1_3_1
conda create --name r-stringr1_3_1 --channel conda-forge --yes r=3.5.1 r-stringr=1.3.1
source activate r-stringr1_3_1
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-stringr1_3_1.yaml

source deactivate

# r-yaml
conda create --name r-yaml --channel conda-forge --yes r=3.5.1 r-yaml=2.2.0
source activate r-yaml
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/r-yaml.yaml

source deactivate

# samtools1_9
conda create --name samtools1_9 --channel bioconda --yes samtools=1.9
source activate samtools1_9
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/samtools1_9.yaml

source deactivate

# singularity
conda create --name singularity --channel conda-forge --channel bioconda --yes singularity=3.0.1
source activate singularity
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/singularity.yaml

source deactivate

#WorkflowError:
#Failed to pull singularity image from docker://continuumio/miniconda3:
#WARNING: pull for Docker Hub is not guaranteed to produce the
#WARNING: same image on repeated pull. Use Singularity Registry
#WARNING: (shub://) to pull exactly equivalent images.
#ERROR: You must install squashfs-tools to build images
#ABORT: Aborting with RETVAL=255
#ERROR: pulling container failed!

conda install --channel conda-forge --yes squashfs-tools=4.3

# sniffles, versions on mac behind that on linux
conda create --name sniffles --channel bioconda --yes sniffles=1.0.11
source activate sniffles
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/sniffles.yaml

source deactivate

wget https://github.com/fritzsedlazeck/Sniffles/archive/master.tar.gz -O output/software/Sniffles.tar.gz
tar xzvf output/software/Sniffles.tar.gz -C output/software
cd output/software/Sniffles-master/
mkdir -p build/
cd build/
cmake ..
make
cd ../bin/sniffles*
./sniffles

# star2_7_0f
conda create --name star2_7_0f --channel conda-forge --channel bioconda --yes star=2.7.0f
source activate star2_7_0f
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/star2_7_0f.yaml

source deactivate

# star2_7_0b
conda create --name star2_7_0b --channel conda-forge --channel bioconda --yes star=2.7.0b
source activate star2_7_0b
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/star2_7_0b.yaml

source deactivate

# subread1_6_3
conda create --name subread1_6_3 --channel conda-forge --channel bioconda --yes subread=1.6.3
source activate subread1_6_3
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/subread1_6_3.yaml

source deactivate

# telseq
conda create --name telseq --channel conda-forge --channel bioconda --yes telseq=0.0.2
source activate telseq
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/telseq.yaml

source deactivate

# ucsc-axtchain
conda create --name ucsc-axtchain --channel conda-forge --channel bioconda --yes ucsc-axtchain=377
source activate ucsc-axtchain
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/ucsc-axtchain.yaml

source deactivate

# ucsc-chainmergesort_ucsc-chainsplit
conda create --name ucsc-chainmergesort_ucsc-chainsplit --channel conda-forge --channel bioconda --yes ucsc-chainmergesort=377 ucsc-chainsplit=377
source activate ucsc-chainmergesort_ucsc-chainsplit
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/ucsc-chainmergesort_ucsc-chainsplit.yaml

source deactivate

# ucsc-chainnet
conda create --name ucsc-chainnet --channel conda-forge --channel bioconda --yes ucsc-chainnet=377
source activate ucsc-chainnet
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/ucsc-chainnet.yaml

source deactivate

# ucsc-fasplit
conda create --name ucsc-fasplit --channel bioconda --yes ucsc-fasplit=377
source activate ucsc-fasplit
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/ucsc-fasplit.yaml

source deactivate

# ucsc-fatotwobit
conda create --name ucsc-fatotwobit --channel bioconda --yes ucsc-fatotwobit=377
source activate ucsc-fatotwobit
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/ucsc-fatotwobit.yaml

source deactivate

# ucsc-liftup
conda create --name ucsc-liftup --channel bioconda --yes ucsc-liftup=377
source activate ucsc-liftup
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/ucsc-liftup.yaml

source deactivate

# ucsc-twobitinfo
conda create --name ucsc-twobitinfo --channel bioconda --yes ucsc-twobitinfo=377
conda activate ucsc-twobitinfo
conda env export --no-builds | sed '$d' | sed '$d' > envs/snakemake/ucsc-twobitinfo.yaml

source deactivate

# wget
conda create --name wget --channel conda-forge --yes wget=1.20.1
source activate wget
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/wget.yaml

source deactivate
