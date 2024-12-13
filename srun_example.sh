srun -n 1 --gpus=1 compute-sanitizer --tool memcheck ./finite_element.cuda 100 0.21 out/gpuattempt_minus/
 srun -n 1 --gpus=1 nsys profile -t cuda --stats=true ./finite_element.cuda 100 0.21 out/gpuattempt_minus/
