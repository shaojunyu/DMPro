cd /home/syu/data1/AutoMethylate

for chr in {1..22}; do
    echo "${chr}"
    srun -c 4 Rscript --vanilla src/CUE/CUE.R ${chr} tmp/AnkeHuels_405K_EPIC res/AnkeHuels_405K_EPIC_baseline > tmp/log/CUE_${chr}.log 2>&1 &
done



for chr in {1..22}; do
    echo "${chr}"
    srun -c 4 Rscript --vanilla src/CUE/CUE.R ${chr} > tmp/log/CUE_${chr}.log 2>&1 &
done