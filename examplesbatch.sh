sbatch --time=00:00:10 --gpus=1 --wrap="./finite_element.cuda 100 0.46 out/gpuattempt_minus/" 
