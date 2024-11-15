mkdir data/XY_MC/GPU_Results

mkdir data/XY_MC/GPU_Results/R_3/

mkdir data/XY_MC/GPU_Results/R_3/long_view_deep

resultFolder=data/XY_MC/GPU_Results/R_3/long_view_deep

rm slurm*.out

for L in 400
do
  for J in 0.15 0.18 0.21 0.24 0.27
  do
  sbatch --time=2-0:0 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
  done
done

for L in 500
do
  for J in 0.12 0.15 0.18 0.21 0.24
  do
  sbatch --time=4-0:0 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
  done
done


for L in 600
do
  for J in 0.12 0.15 0.18 0.21
  do
  sbatch --time=05:10:00 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
  done
done