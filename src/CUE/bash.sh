cd /home/syu/data1/AutoMethylate
Rscript --vanilla src/CUE/baseline.R 1 tmp/AnkeHuels_405K_EPIC res/AnkeHuels_405K_EPIC_baseline

for chr in {1..22}; do
    echo "${chr}"
    srun -c 4 Rscript --vanilla src/CUE/CUE.R ${chr} tmp/AnkeHuels_405K_EPIC res/AnkeHuels_405K_EPIC_baseline > tmp/log/CUE_${chr}.log 2>&1 &
done

# Anke
for chr in {1..22}; do
    echo "${chr}"
    srun -c 8 Rscript --vanilla src/CUE/baseline.R ${chr} tmp/AnkeHuels_405K_EPIC res/AnkeHuels_405K_EPIC_baseline > tmp/log/Anke_baseline_${chr}.log 2>&1 &
done

# Alica
Rscript --vanilla src/CUE/baseline.R 20 tmp/AliciaKSmith_450K_EPIC res/AliciaKSmith_450K_EPIC_baseline
for chr in {1..22}; do
    echo "${chr}"
    srun -c 8 Rscript --vanilla src/CUE/baseline.R ${chr} tmp/AliciaKSmith_450K_EPIC res/AliciaKSmith_450K_EPIC_baseline > tmp/log/Alica_baseline_${chr}.log 2>&1 &
done

# CUE
for chr in {1..22}; do
    echo "${chr}"
    srun -c 4 Rscript --vanilla src/CUE/CUE.R ${chr} tmp/AliciaKSmith_450K_EPIC res/AliciaKSmith_450K_EPIC_baseline
done

for chr in {1..22}; do
    echo "${chr}"
    srun -c 4 Rscript --vanilla src/CUE/CUE.R ${chr} tmp/AnkeHuels_405K_EPIC res/AnkeHuels_405K_EPIC_baseline &
done