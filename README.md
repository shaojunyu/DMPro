# DMPro
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
- 171 samples from Anke Huels (anke.huels@emory.edu) (405K and EPIC)
- 123 samples from Alicia K. Smith (aksmit3@emory.edu) (450K and EPIC)
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

## model training, Anke Huels's data
conda activate automethy
for chr in {1..22}; do
    echo $chr
    python src/AE.py \
        --mode train \
        --platform EPIC \
        --chr $chr \
        --data_prefix tmp/AnkeHuels_405K_EPIC \
        --model_path res/AnkeHuels_405K_EPIC_models \
        --bin_size_file src/Anke_bin_size.csv > tmp/log/AnkeHuels_train_EPIC_chr$chr.log 2>&1
done
for chr in {1..22}; do
    echo $chr
    python src/AE.py \
        --mode train \
        --platform HM450 \
        --chr $chr \
        --data_prefix tmp/AnkeHuels_405K_EPIC \
        --model_path res/AnkeHuels_405K_EPIC_models \
        --bin_size_file src/Anke_bin_size.csv > tmp/log/AnkeHuels_train_HM450_chr$chr.log 2>&1
done

## model testing, Anke Huels's data
for chr in {1..22}; do
    python src/AE.py \
        --mode test \
        --platform EPIC \
        --chr $chr \
        --data_prefix tmp/AnkeHuels_405K_EPIC \
        --model_path res/AnkeHuels_405K_EPIC_models \
        --bin_size_file src/Anke_bin_size.csv
done

for chr in {1..22}; do
    python src/AE.py \
        --mode test \
        --platform HM450 \
        --chr $chr \
        --data_prefix tmp/AnkeHuels_405K_EPIC \
        --model_path res/AnkeHuels_405K_EPIC_models \
        --bin_size_file src/Anke_bin_size.csv
done

# model training for Alicia Smith's data
for chr in {1..22}; do
    echo $chr
    python src/AE.py \
        --mode train \
        --platform EPIC \
        --chr $chr \
        --data_prefix tmp/AliciaKSmith_450K_EPIC \
        --model_path res/AliciaKSmith_450K_EPIC_models \
        --bin_size_file src/Alicia_bin_size.csv > tmp/log/Alicia_train_EPIC_chr$chr.log 2>&1
done
for chr in {1..22}; do
    echo $chr
    python src/AE.py \
        --mode train \
        --platform HM450 \
        --chr $chr \
        --data_prefix tmp/AliciaKSmith_450K_EPIC \
        --model_path res/AliciaKSmith_450K_EPIC_models \
        --bin_size_file src/Alicia_bin_size.csv > tmp/log/Alicia_train_HM450_chr$chr.log 2>&1
done

# Testing
for chr in {1..22}; do
    python src/AE.py \
        --mode test \
        --platform EPIC \
        --chr $chr \
        --data_prefix tmp/AliciaKSmith_450K_EPIC \
        --model_path res/AliciaKSmith_450K_EPIC_models \
        --bin_size_file src/Alicia_bin_size.csv
done
for chr in {1..22}; do
    python src/AE.py \
        --mode test \
        --platform HM450 \
        --chr $chr \
        --data_prefix tmp/AliciaKSmith_450K_EPIC \
        --model_path res/AliciaKSmith_450K_EPIC_models \
        --bin_size_file src/Alicia_bin_size.csv
done

## Encode Alicia Smith's data
for chr in {1..22}; do
    python src/AE.py \
        --mode encode \
        --platform HM450 \
        --chr $chr \
        --data_prefix tmp/AliciaKSmith_450K_EPIC \
        --model_path res/AliciaKSmith_450K_EPIC_models \
        --bin_size_file src/Alicia_bin_size.csv \
        --output_path res/AliciaKSmith_450K_EPIC_encoded
done
for chr in {1..22}; do
    python src/AE.py \
        --mode encode \
        --platform EPIC \
        --chr $chr \
        --data_prefix tmp/AliciaKSmith_450K_EPIC \
        --model_path res/AliciaKSmith_450K_EPIC_models \
        --bin_size_file src/Alicia_bin_size.csv \
        --output_path res/AliciaKSmith_450K_EPIC_encoded
done

## Encode Anke Huels's data
for chr in {1..22}; do
    python src/AE.py \
        --mode encode \
        --platform EPIC \
        --chr $chr \
        --data_prefix tmp/AnkeHuels_405K_EPIC \
        --model_path res/AnkeHuels_405K_EPIC_models \
        --bin_size_file src/Anke_bin_size.csv \
        --output_path res/AnkeHuels_405K_EPIC_encoded
done
for chr in {1..22}; do
    python src/AE.py \
        --mode encode \
        --platform HM450 \
        --chr $chr \
        --data_prefix tmp/AnkeHuels_405K_EPIC \
        --model_path res/AnkeHuels_405K_EPIC_models \
        --bin_size_file src/Anke_bin_size.csv \
        --output_path res/AnkeHuels_405K_EPIC_encoded
done


### NMF
for chr in {1..22}; do
    echo $chr
    python src/NMF.py --chr $chr \
        --encoding_path res/AnkeHuels_405K_EPIC_encoded \
        --data_path tmp/AnkeHuels_405K_EPIC --out_path res/AnkeHuels_405K_EPIC_NMF
done

for chr in {1..1}; do
    echo $chr
    python src/NMF.py --chr $chr \
        --encoding_path res/AliciaKSmith_450K_EPIC_encoded \
        --data_path tmp/AliciaKSmith_450K_EPIC --out_path res/AliciaKSmith_450K_EPIC_NMF
done

### decode
export chr=1
python src/AE.py --mode decode --platform EPIC --chr $chr --model_path res/AnkeHuels_405K_EPIC_models --encoding_file res/AnkeHuels_405K_EPIC_NMF/chr${chr}_reconstruct_new.csv --data_prefix tmp/AnkeHuels_405K_EPIC --bin_size_file src/Anke_bin_size.csv --decode_file res/AnkeHuels_405K_EPIC_decoded/chr${chr}_decoded_new.csv

for chr in {1..22}; do
    echo $chr
    python src/AE.py --mode decode --platform EPIC --chr $chr --model_path res/AnkeHuels_405K_EPIC_models --encoding_file res/AnkeHuels_405K_EPIC_NMF/chr${chr}_reconstruct_new.csv --data_prefix tmp/AnkeHuels_405K_EPIC --bin_size_file src/Anke_bin_size.csv --decode_file res/AnkeHuels_405K_EPIC_decoded/chr${chr}_decoded_new.csv
done

for chr in {1..22}; do
    echo $chr
    python src/AE.py --mode decode --platform EPIC \
    --chr $chr --model_path res/AliciaKSmith_450K_EPIC_models \
    --encoding_file res/AliciaKSmith_450K_EPIC_NMF/chr${chr}_reconstruct.csv \
    --data_prefix tmp/AliciaKSmith_450K_EPIC --bin_size_file src/Alicia_bin_size.csv \
    --decode_file res/AliciaKSmith_450K_EPIC_decoded/chr${chr}_decoded.csv
done

# transfer learning

# encode Alicia Smith's data using Anke Huels's model
export chr=20

# test model
python src/AE.py --mode test --platform EPIC --chr $chr --model_path res/AnkeHuels_405K_EPIC_models --data_prefix tmp/AliciaKSmith_450K_EPIC --bin_size_file src/Alicia_bin_size.csv

python src/AE.py --mode test --platform EPIC --chr $chr --model_path res/AliciaKSmith_450K_EPIC_models --data_prefix tmp/AliciaKSmith_450K_EPIC --bin_size_file src/Alicia_bin_size.csv


python src/AE.py --mode encode --platform EPIC --chr $chr --model_path res/AnkeHuels_405K_EPIC_models --data_prefix tmp/AliciaKSmith_450K_EPIC --bin_size_file src/Alicia_bin_size.csv --output_path res/AliciaKSmith_450K_EPIC_encoded_AH
python src/AE.py --mode encode --platform HM450 --chr $chr --model_path res/AnkeHuels_405K_EPIC_models --data_prefix tmp/AliciaKSmith_450K_EPIC --bin_size_file src/Alicia_bin_size.csv --output_path res/AliciaKSmith_450K_EPIC_encoded_AH

# NMF
python src/NMF.py --chr $chr --encoding_path res/AliciaKSmith_450K_EPIC_encoded_AH --data_path tmp/AliciaKSmith_450K_EPIC --out_path res/AliciaKSmith_450K_EPIC_NMF_AH

# decode
python src/AE.py --mode decode --platform EPIC --chr $chr --model_path res/AnkeHuels_405K_EPIC_models --encoding_file res/AliciaKSmith_450K_EPIC_NMF_AH/chr${chr}_reconstruct.csv --data_prefix tmp/AliciaKSmith_450K_EPIC --bin_size_file src/Alicia_bin_size.csv --decode_file res/AliciaKSmith_450K_EPIC_decoded_AH/chr${chr}_decoded.csv
```
