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
- EPIC Series: https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&platform=21145&display=20&zsort=samples

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

## Usage
```bash

# train 450K
python src/AE.py --mode train --platform HM450 --bin_size 4000  --chr 20 --epochs 200 --step_size 1000
python src/AE.py --mode test --platform HM450 --bin_size 4000  --chr 20 --step_size 1000
python src/AE.py --mode encode --platform HM450 --bin_size 4000  --chr 20


# train EPIC
python src/AE.py --mode train --platform EPIC --bin_size 8000 --chr 20 --epochs 200 --step_size 1000
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 20 --step_size 1000
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --data_prefix tmp/processed/EPIC --chr 20
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --data_prefix tmp/processed/EPIC --chr 20


for chr in {1..22}; do
    export CUDA_VISIBLE_DEVICES=0
    echo $chr
    python src/AE.py --mode train --platform HM450 --bin_size 4000  --chr $chr --epochs 200 --step_size 100 > tmp/log/train_HM450_chr$chr.log 2>&1
done
for chr in {1..22}; do
    export CUDA_VISIBLE_DEVICES=1
    echo $chr
    python src/AE.py --mode train --platform EPIC --bin_size 8000  --chr $chr --epochs 300 --step_size 100 > tmp/log/train_EPIC_chr$chr.log 2>&1
done
export chr=20;
python src/AE.py --mode train --platform HM450 --bin_size 4000  --chr $chr --epochs 200 --step_size 100 > tmp/log/train_HM450_chr$chr.log 2>&1
export chr=20;
python src/AE.py --mode train --platform EPIC --bin_size 8000  --chr $chr --epochs 200 --step_size 100 > tmp/log/train_EPIC_chr$chr.log 2>&1


#chr1
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 1  --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 1  --step_size 3200
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 1  --step_size 3200
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 1  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 1  --step_size 5000


#chr2
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 2 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 2  --step_size 3000
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 2  --step_size 3000
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 2  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 2  --step_size 5000

#chr3
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 3 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 3  --step_size 3000
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 3  --step_size 3000
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 3  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 3  --step_size 5000

#chr4
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 4 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 4  --step_size 3200
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 4  --step_size 3200
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 4  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 4  --step_size 5000

#chr5
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 5 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 5  --step_size 3000
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 5  --step_size 3000
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 5  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 5  --step_size 5000

#chr6
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 6 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 6  --step_size 3800
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 6  --step_size 3800
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 6  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 6  --step_size 5000

#chr7
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 7 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 7  --step_size 3600
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 7  --step_size 3600
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 7  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 7  --step_size 5000

#chr8
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 8 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 8  --step_size 3200
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 8  --step_size 3200
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 8  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 8  --step_size 5000

#chr9
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 9 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 9  --step_size 2600
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 9  --step_size 2600
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 9  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 9  --step_size 5000

#chr10
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 10 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 10  --step_size 3200
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 10  --step_size 3200
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 10  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 10  --step_size 5000

#chr11
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 11 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 11  --step_size 3400
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 11  --step_size 3400
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 11  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 11  --step_size 5000

#chr12
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 12 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 12  --step_size 3000
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 12  --step_size 3000
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 12  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 12  --step_size 5000


#chr13
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 13 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 13  --step_size 3000
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 13  --step_size 3000
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 13  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 13  --step_size 5000

#chr14
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 14 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 14  --step_size 3000
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 14  --step_size 3000
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 14  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 14  --step_size 5000

#chr15
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 15 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 15  --step_size 3000
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 15  --step_size 3000
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 15  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 15  --step_size 5000

#chr16
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 16 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 16  --step_size 3000
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 16  --step_size 3000
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 16  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 16  --step_size 5000

#chr17
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 17 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 17  --step_size 3200
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 17  --step_size 3200
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 17  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 17  --step_size 5000

#chr18
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 18 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 18  --step_size 1000
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 18  --step_size 1000
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 18  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 18  --step_size 5000

#chr19
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 19 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 19  --step_size 3400
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 19  --step_size 3400
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 19  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 19  --step_size 5000

#chr20
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 20 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 20  --step_size 2800
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 20  --step_size 2800
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 20  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 20  --step_size 5000


#chr21
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 21 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 21  --step_size 2400
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 21  --step_size 2400
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 21  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 21  --step_size 5000

#chr22
python src/AE.py --mode test --platform EPIC --bin_size 8000 --chr 22 --step_size 5000
python src/AE.py --mode test --platform HM450 --bin_size 4000 --chr 22  --step_size 2400
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 22  --step_size 2400
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --chr 22  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --chr 22  --step_size 5000



# chr18 450k
python src/AE.py --mode train --platform 450k --bin_size 4000 --data_prefix tmp/NMF/HM450 --chr 18 --epochs 200
python src/AE.py --mode test --platform 450k --bin_size 4000 --data_prefix tmp/NMF/HM450 --chr 18 --step_size 1000
python src/AE.py --mode encode --platform HM450 --bin_size 4000 --chr 20  --step_size 1000

# chr18 EPIC
python src/AE.py --mode train --platform EPIC --bin_size 8000 --data_prefix tmp/NMF/EPIC --chr 18 --epochs 300
python src/AE.py --mode test --platform EPIC --bin_size 8000 --data_prefix tmp/NMF/EPIC --chr 18 --step_size 5000
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --data_prefix tmp/NMF/EPIC --chr 18  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --data_prefix tmp/NMF/EPIC --chr 18  --step_size 5000

# chr19 450k
python src/AE.py --mode train --platform 450k --bin_size 4000 --data_prefix tmp/NMF/HM450 --chr 19 --step_size 1000 --epochs 200
python src/AE.py --mode test --platform 450k --bin_size 4000 --data_prefix tmp/NMF/HM450 --chr 19 --step_size 4000
python src/AE.py --mode encode --platform 450k --bin_size 4000 --data_prefix tmp/NMF/HM450 --chr 19  --step_size 4000

# chr19 EPIC
python src/AE.py --mode train --platform EPIC --bin_size 8000 --data_prefix tmp/NMF/EPIC --chr 19 --step_size 2000 --epochs 300
python src/AE.py --mode test --platform EPIC --bin_size 8000 --data_prefix tmp/NMF/EPIC --chr 19 --step_size 5000
python src/AE.py --mode encode --platform EPIC --bin_size 8000 --data_prefix tmp/NMF/EPIC --chr 19  --step_size 5000
python src/AE.py --mode decode --platform EPIC --bin_size 8000 --data_prefix tmp/NMF/EPIC --chr 19  --step_size 5000

for chr in {1..22}; do
    srun -c 10 Rscript src/CUE/CUE.R $chr &
done

for chr in {1..22}; do
    srun -c 2 Rscript src/CUE/CUE.R $chr &
done

```