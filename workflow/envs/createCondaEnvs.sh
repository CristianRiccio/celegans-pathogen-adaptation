# Create conda environments

# nick
conda create --name nick --channel conda-forge --channel bioconda --yes git=2.21.0 r=3.5.1 snakemake=5.4.0
source activate nick
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/nick.yaml

source deactivate

# blat
conda create --name blat --channel bioconda --yes blat=36
source activate blat
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/blat.yaml

source deactivate

# bcftools
conda create --name bcftools --channel conda-forge --channel bioconda --yes bcftools=1.9
source activate bcftools
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/bcftools.yaml

source deactivate

# bedtools
conda create --name bedtools --channel conda-forge --channel bioconda --yes bedtools=2.27.1
source activate bedtools
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/bedtools.yaml

source deactivate

bioconductor-biomart_bioconductor-deseq2_r-stringr_r-tidyr
conda create --name bioconductor-biomart_bioconductor-deseq2_r-stringr_r-tidyr --channel conda-forge --yes bioconductor-biomart=2.38.0 bioconductor-deseq2=1.22.1 r-stringr=1.3.1 r-tidyr=0.8.3
source activate bioconductor-biomart_bioconductor-deseq2_r-stringr_r-tidyr
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/bioconductor-biomart_bioconductor-deseq2_r-stringr_r-tidyr.yaml

source deactivate

# bioconductor-deseq2
conda create --name bioconductor-deseq2 --channel conda-forge --channel bioconda --yes bioconductor-deseq2=1.22.1
source activate bioconductor-deseq2
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/bioconductor-deseq2.yaml

source deactivate

# r=3.6_bioconductor-deseq2=1.26.0_r-stringr=1.4.0
conda create --name r=3.6_bioconductor-deseq2=1.26.0_r-stringr=1.4.0 --channel conda-forge --channel bioconda --channel r --yes r=3.6 bioconductor-deseq2=1.26.0 r-stringr=1.4.0
source activate r=3.6_bioconductor-deseq2=1.26.0_r-stringr=1.4.0
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r=3.6_bioconductor-deseq2=1.26.0_r-stringr=1.4.0.yaml

source deactivate

# bioconductor-deseq2_r-dplyr_r-stringr
conda create --name bioconductor-deseq2_r-dplyr_r-stringr --channel conda-forge --channel bioconda --yes bioconductor-deseq2=1.22.1 r=3.5.1 r-dplyr=0.8.0 r-stringr=1.3.1
source activate bioconductor-deseq2_r-dplyr_r-stringr
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/bioconductor-deseq2_r-dplyr_r-stringr.yaml

source deactivate

# bioconductor-deseq2_r-plotly_r-stringr
conda create --name bioconductor-deseq2_r-plotly_r-stringr --channel conda-forge --channel bioconda --yes bioconductor-deseq2=1.22.1 r=3.5.1 r-plotly=4.8.0 r-stringr=1.3.1
source activate bioconductor-deseq2_r-plotly_r-stringr
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/bioconductor-deseq2_r-plotly_r-stringr.yaml

source deactivate

# bioconductor-deseq2_r-dplyr_r-pheatmap_r-stringr
conda create --name bioconductor-deseq2_r-dplyr_r-pheatmap_r-stringr --channel conda-forge --channel bioconda --yes bioconductor-deseq2=1.22.1 r=3.5.1 r-dplyr=0.8.0 r-pheatmap=1.0.12 r-stringr=1.3.1
source activate bioconductor-deseq2_r-dplyr_r-pheatmap_r-stringr
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/bioconductor-deseq2_r-dplyr_r-pheatmap_r-stringr.yaml

source deactivate

# bioconductor-rtracklayer
conda create --name bioconductor-rtracklayer --channel conda-forge --channel bioconda --channel r --yes bioconductor-rtracklayer=1.42.1 r=3.5.1
source activate bioconductor-rtracklayer
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/bioconductor-rtracklayer.yaml

source deactivate

# bioconductor-rsamtools_r-seqinr_r-tidyr
conda create --name bioconductor-rsamtools_r-seqinr_r-tidyr --channel conda-forge --channel bioconda --yes bioconductor-rsamtools=1.34.0 r=3.5.1 r-seqinr=3.4_5 r-tidyr=0.8.2
conda activate bioconductor-rsamtools_r-seqinr_r-tidyr
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/bioconductor-rsamtools_r-seqinr_r-tidyr.yaml

source deactivate

# bioconductor-rtracklayer_r-dplyr_r-ggplot2
conda create --name bioconductor-rtracklayer_r-dplyr_r-ggplot2 --channel conda-forge --channel bioconda --yes r=3.5.1 bioconductor-rtracklayer=1.42.1 r-dplyr=0.8.0 r-ggplot2=3.1.0
source activate bioconductor-rtracklayer_r-dplyr_r-ggplot2
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/bioconductor-rtracklayer_r-dplyr_r-ggplot2.yaml

source deactivate

# bwa
conda create --name bwa --channel bioconda --yes bwa=0.7.17
source activate bwa
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/bwa.yaml

source deactivate

# canu, version on Linux ahead of Mac
conda create --name canu --channel bioconda --yes canu=1.8
source activate canu
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/canu.yaml

source deactivate

# csvkit, openpyxl 2.6.0 breaks csvkit
conda create --name csvkit --channel conda-forge --yes csvkit=1.0.3 openpyxl=2.5.0
source activate csvkit
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/csvkit.yaml

source deactivate

# cufflinks2_2_1
conda create --name cufflinks2_2_1 --channel conda-forge --channel bioconda --channel r --yes cufflinks=2.2.1
source activate cufflinks2_2_1
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/cufflinks2_2_1.yaml

source deactivate

# cutadapt1_18
conda create --name cutadapt1_18 --channel bioconda --yes cutadapt=1.18
source activate cutadapt1_18
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/cutadapt1_18.yaml

source deactivate

# fastqc0_11_8
conda create --name fastqc0_11_8 --channel conda-forge --channel bioconda --yes fastqc=0.11.8
source activate fastqc0_11_8
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/fastqc0_11_8.yaml

source deactivate

# fastx_toolkit
conda create --name fastx_toolkit --channel bioconda --yes fastx_toolkit=0.0.14
source activate fastx_toolkit
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/fastx_toolkit.yaml

source deactivate

# hisat2v2_1_0
conda create --name hisat2v2_1_0 --channel bioconda --channel conda-forge --yes hisat2=2.1.0
source activate hisat2v2_1_0
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/hisat2v2_1_0.yaml

source deactivate

# htseq
conda create --name htseq --channel conda-forge --channel bioconda --yes htseq=0.11.2
conda activate htseq
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/htseq.yaml

source deactivate


# kallisto0_45_1
conda create --name kallisto0_45_1 --channel conda-forge --channel bioconda --yes kallisto=0.45.1
source activate kallisto0_45_1
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/kallisto0_45_1.yaml

source deactivate

# minimap2
conda create --name minimap2 --channel bioconda --yes minimap2=2.15
source activate minimap2
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/minimap2.yaml

source deactivate

# multiqc1_7
conda create --name multiqc1_7 --channel conda-forge --channel bioconda --yes multiqc=1.7
source activate multiqc1_7
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/multiqc1_7.yaml

source deactivate

# mummer=3.23
conda create --name mummer=3.23 --channel bioconda --yes mummer=3.23
conda activate mummer=3.23
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/mummer=3.23.yaml

conda deactivate

# picard
conda create --name picard --channel bioconda --yes picard=2.18.26
source activate picard
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/picard.yaml

source deactivate

# python3_7_3
conda create --name python3_7_3 --channel conda-forge --yes python=3.7.3
source activate python3_7_3
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/python3_7_3.yaml

source deactivate

# plotQC_N2
conda create --name plotQC_N2 --channel conda-forge --yes r r-dplyr r-ggplot2
source activate plotQC_N2
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/plotQC_N2.yaml

source deactivate

# r3_5_1
conda create --name r --channel conda-forge --yes r=3.5.1
source activate r3_5_1
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r3_5_1.yaml

# r-all
conda create --name r-all --channel conda-forge --yes  bioconductor-deseq2=1.22.1 bioconductor-rtracklayer=1.42.1 bioconductor-rsamtools=1.34.0 r=3.5.1 r-dplyr=0.8.1 r-ggplot2=3.1.0 r-plotly=4.8.0 r-seqinr=3.4_5 r-stringr=1.4.0 r-tidyr=0.8.3 r-yaml=2.2.0
source activate r-all
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r-all.yaml

source deactivate



# r=3.6_bioconductor-all=1.28.0
conda create --name r=3.6_bioconductor-all=1.28.0 --channel conda-forge --channel bioconda --yes r=3.6 bioconductor-all=1.28.0
source activate r=3.6_bioconductor-all=1.28.0
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r=3.6_bioconductor-all=1.28.0.yaml

source deactivate

# r=3.6_bioconductor-annotationdbi=1.48.0
conda create --name r=3.6_bioconductor-annotationdbi=1.48.0 --channel conda-forge --channel bioconda --yes r=3.6 bioconductor-annotationdbi=1.48.0
source activate r=3.6_bioconductor-annotationdbi=1.48.0
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r=3.6_bioconductor-annotationdbi=1.48.0.yaml

source deactivate

# r=3.6_bioconductor-all=1.28.0_bioconductor-annotationdbi=1.48.0_bioconductor-go.db=3.10.0_bioconductor-rgraphviz=2.30.0_bioconductor-topgo=2.37.0
conda create --name r=3.6_bioconductor-all=1.28.0_bioconductor-annotationdbi=1.48.0_bioconductor-go.db=3.10.0_bioconductor-rgraphviz=2.30.0_bioconductor-topgo=2.37.0 --channel conda-forge --channel bioconda --yes r=3.6 bioconductor-all bioconductor-annotationdbi=1.48.0 bioconductor-go.db=3.10.0 bioconductor-rgraphviz=2.30.0 bioconductor-topgo=2.37.0
source activate r=3.6_bioconductor-all=1.28.0_bioconductor-annotationdbi=1.48.0_bioconductor-go.db=3.10.0_bioconductor-rgraphviz=2.30.0_bioconductor-topgo=2.37.0
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r=3.6_bioconductor-all=1.28.0_bioconductor-annotationdbi=1.48.0_bioconductor-go.db=3.10.0_bioconductor-rgraphviz=2.30.0_bioconductor-topgo=2.37.0.yaml

source deactivate

# r=3.6_bioconductor-biomart=2.42.0_r-tidyr=1.0.2
conda create --name r=3.6_bioconductor-biomart=2.42.0_r-tidyr=1.0.2 --channel conda-forge --yes r=3.6 bioconductor-biomart=2.42.0 r-tidyr=1.0.2
source activate r=3.6_bioconductor-biomart=2.42.0_r-tidyr=1.0.2
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r=3.6_bioconductor-biomart=2.42.0_r-tidyr=1.0.2.yaml

source deactivate



# r=3.6_bioconductor-go.db=3.10.0
conda create --name r=3.6_bioconductor-go.db=3.10.0 --channel conda-forge --channel bioconda --yes r=3.6 bioconductor-go.db=3.10.0
source activate r=3.6_bioconductor-go.db=3.10.0
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r=3.6_bioconductor-go.db=3.10.0.yaml

source deactivate

# r=3.6_bioconductor-topgo=2.37.0
conda create --name r=3.6_bioconductor-topgo=2.37.0 --channel conda-forge --channel bioconda --yes r=3.6 bioconductor-topgo=2.37.0
source activate r=3.6_bioconductor-topgo=2.37.0
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r=3.6_bioconductor-topgo=2.37.0.yaml

source deactivate

# r-dplyr
conda create --name r-dplyr --channel conda-forge --yes r=3.5.1 r-dplyr=0.8.1
source activate r-dplyr
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r-dplyr.yaml

source deactivate

# r-dplyr_r-stringr
conda create --name r-dplyr_r-stringr --channel conda-forge --yes r=3.5.1 r-dplyr=0.8.1  r-stringr=1.4.0
source activate r-dplyr_r-stringr
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r-dplyr_r-stringr.yaml

source deactivate

# r-ggplot2
conda create --name r-ggplot2 --channel conda-forge --yes r=3.5.1 r-ggplot2=3.1.0
source activate r-ggplot2
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r-ggplot2.yaml

source deactivate

# r-ggplot2_r-plotly_r-stringr_r-tidyr
conda create --name r-ggplot2_r-plotly_r-stringr_r-tidyr --channel conda-forge --yes r=3.5.1 r-ggplot2=3.1. r-plotly=4.8.0 r-stringr=1.3.1 r-tidyr=0.8.3
source activate r-ggplot2_r-plotly_r-stringr_r-tidyr
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r-ggplot2_r-plotly_r-stringr_r-tidyr.yaml

source deactivate

# r-plotly_r-tidyverse
conda create --name r-ggplot2 --channel conda-forge --yes r=3.5.1 r-plotly=4.8.0 r-tidyverse=1.2.1
source activate r-plotly_r-tidyverse
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r-plotly_r-tidyverse.yaml

source deactivate

# r-seqinr
conda create --name r-seqinr --channel conda-forge --yes r=3.5.1 r-seqinr=3.4_5
source activate r-seqinr
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r-seqinr.yaml

source deactivate

# r-stringr1_3_1
conda create --name r-stringr1_3_1 --channel conda-forge --yes r=3.5.1 r-stringr=1.3.1
source activate r-stringr1_3_1
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r-stringr1_3_1.yaml

source deactivate

# r-yaml
conda create --name r-yaml --channel conda-forge --yes r=3.5.1 r-yaml=2.2.0
source activate r-yaml
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/r-yaml.yaml

source deactivate

# samtools1_9
conda create --name samtools1_9 --channel bioconda --yes samtools=1.9
source activate samtools1_9
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/samtools1_9.yaml

source deactivate

# conda install --channel conda-forge --yes squashfs-tools=4.3

# star2_7_0f
conda create --name star2_7_0f --channel conda-forge --channel bioconda --yes star=2.7.0f
source activate star2_7_0f
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/star2_7_0f.yaml

source deactivate

# star2_7_0b
conda create --name star2_7_0b --channel conda-forge --channel bioconda --yes star=2.7.0b
source activate star2_7_0b
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/star2_7_0b.yaml

source deactivate

# subread1_6_3
conda create --name subread1_6_3 --channel conda-forge --channel bioconda --yes subread=1.6.3
source activate subread1_6_3
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/subread1_6_3.yaml

source deactivate

# subread=1.6.3
conda create --name subread1.6.3 --channel conda-forge --channel bioconda --yes subread=1.6.3
source activate subread1.6.3
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/subread1.6.3.yaml

source deactivate

# telseq
conda create --name telseq --channel conda-forge --channel bioconda --yes telseq=0.0.2
source activate telseq
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/telseq.yaml

source deactivate

# wget
conda create --name wget --channel conda-forge --yes wget=1.20.1
source activate wget
conda env export --no-builds | sed '$d' | sed '$d' > workflow/envs/conda/wget.yaml

source deactivate
