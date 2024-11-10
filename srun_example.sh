srun -n 1 --gpus=1 compute-sanitizer --tool memcheck ./finite_element.cuda 100 0.21 out/gpuattempt_minus/
