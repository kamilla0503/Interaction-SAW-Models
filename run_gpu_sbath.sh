mkdir data/XY_MC/GPU_Results

mkdir data/XY_MC/GPU_Results/R_3/

mkdir data/XY_MC/GPU_Results/R_3/long_view

resultFolder=data/XY_MC/GPU_Results/R_3/long_view


for L in 100 200
do
  for J in 0.25 0.35
  do
  sbatch --time=01:40:00 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
  done
done


for L in 300 400
do
  for J in 0.2 0.3 0.4 0.25 0.35
  do
  sbatch --time=02:10:00 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
  done
done