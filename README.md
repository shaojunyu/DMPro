# AutoMethylate
DNA methylation imputation tool powered by the fusion of autoencoder and NMF algorithms

## Information about the methylation array
| Platform | Short Name| GPL Number | Number of CpG sites | Release Date | URL |
| --- | --- | --- | --- | --- | --- |
Illumina HumanMethylation27 BeadChip | 27K | GPL8490 | 27,578 | Apr 27, 2009 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8490
Illumina HumanMethylation450 BeadChip | 450K | GPL13534 | 485,577 | May 13, 2011 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534
Illumina Infinium HumanMethylation850 BeadChip | 850K | GPL23976 | 866,895 | Sep 01, 2017 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL23976
Illumina Infinium MethylationEPIC | EPIC | GPL21145 | 868,564 | Nov 16, 2015 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21145
Illumina Infinium MethylationEPIC version 2 | EPICv2 | GPL33022 | 937,691 | Jan 15, 2023 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL33022

## Data
- 171 samples
- data are orginized by chromosomes

## Data for training autoencoder
- https://www.ncbi.nlm.nih.gov/geo/browse/?view=platforms
- 450K Series: https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&platform=13534&display=20&zsort=samples

## Runtime Environment
- Python 3.7.3
    - numpy 1.16.4
- R 4.2.0
    - data.table
    - tidyverse

```bash
conda create -n automethy python=3.10
conda activate automethy
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia

```