mkdir data/XY_MC/GPU_Results

mkdir data/XY_MC/GPU_Results/R_3/

mkdir data/XY_MC/GPU_Results/R_3/long_view_short

resultFolder=data/XY_MC/GPU_Results/R_3/long_view_short

#rm slurm*.out

for L in 200
do
  for J in 
  do
  sbatch --time=2-0:0 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
  done
done


# for L in 400
# do
#   for J in 0.246 0.252 0.258 0.264
#   do
#   sbatch --time=2-0:0 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
#   done
# done


# for L in 500
# do
#   for J in 0.246 0.252 0.258 0.264
#   do
#   sbatch --time=4-0:0 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
#   done
# done


# for L in 600
# do
#   for J in 0.246 0.252 0.258 0.264
#   do
#   sbatch --time=6-0:0 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
#   done
# done
