mkdir data/XY_MC/GPU_Results

mkdir data/XY_MC/GPU_Results/R_3/

mkdir data/XY_MC/GPU_Results/R_3/long_view_short

resultFolder=data/XY_MC/GPU_Results/R_3/long_view_short

#rm slurm*.out


for L in 400
do
  for J in 0.277 0.279
  do
  sbatch --time=2-0:0 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
  done
done


for L in 500
do
  for J in 0.268, 0.272, 0.276, 0.28 
  do
  sbatch --time=4-0:0 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
  done
done



# for L in 400
# do
#   for J in 0.28 0.275
#   do
#   sbatch --time=2-0:0 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
#   done
# done

# for L in 300
# do
#   for J in 0.22 0.24 0.26 0.28 0.3
#   do
#   sbatch --time=0-5:0 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
#   done
# done


# for L in 200
# do
#   for J in 0.31 0.305
#   do
#   sbatch --time=0-2:0 --gpus=1 --wrap="./finite_element.cuda $L $J $resultFolder"
#   done
# done

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
