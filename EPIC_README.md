
## Usage
```bash

## model training, EPIC V1 data
conda activate automethy
for chr in {1..22}; do
    echo $chr
    python src/EPIC_AE.py \
        --mode train \
        --platform EPIC_V1 \
        --chr $chr \
        --data_prefix tmp/EPIC_V1_V2 \
        --model_path res/EPIC_V1_V2_models \
        --bin_size_file src/EPIC_V1_V2_bin_size.csv > tmp/log/EPIC_V1_train_EPIC_chr$chr.log 2>&1
done
# EPIC V2
for chr in {1..22}; do
    echo $chr
    python src/EPIC_AE.py \
        --mode train \
        --platform EPIC_V2 \
        --chr $chr \
        --data_prefix tmp/EPIC_V1_V2 \
        --model_path res/EPIC_V1_V2_models \
        --bin_size_file src/EPIC_V1_V2_bin_size.csv > tmp/log/EPIC_V2_train_EPIC_chr$chr.log 2>&1
done


# Testing
for chr in {1..22}; do
    python src/EPIC_AE.py \
        --mode test \
        --platform EPIC_V2 \
        --chr $chr \
        --data_prefix tmp/EPIC_V1_V2 \
        --model_path res/EPIC_V1_V2_models \
        --bin_size_file src/EPIC_V1_V2_bin_size.csv
done



## Encode EPIC V1 data
for chr in {1..22}; do
    python src/EPIC_AE.py \
        --mode encode \
        --platform EPIC_V1 \
        --chr $chr \
        --data_prefix tmp/EPIC_V1_V2 \
        --model_path res/EPIC_V1_V2_models \
        --bin_size_file src/EPIC_V1_V2_bin_size.csv \
        --output_path res/EPIC_V1_V2_encoded &
done
for chr in {1..22}; do
    python src/EPIC_AE.py \
        --mode encode \
        --platform EPIC_V2 \
        --chr $chr \
        --data_prefix tmp/EPIC_V1_V2 \
        --model_path res/EPIC_V1_V2_models \
        --bin_size_file src/EPIC_V1_V2_bin_size.csv \
        --output_path res/EPIC_V1_V2_encoded &
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
export chr=20
python src/EPIC_AE.py --mode decode --platform EPIC_V2 --chr $chr --model_path res/EPIC_V1_V2_models --encoding_file res/EPIC_V1_V2_NMF/chr${chr}_reconstruct_new.csv --data_prefix tmp/EPIC_V1_V2 --bin_size_file src/EPIC_V1_V2_bin_size.csv --decode_file res/EPIC_V1_V2_decoded/chr${chr}_decoded_new.csv

```