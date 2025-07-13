module load CUDA
module load cmake/3.21.3
rm -rf build; mkdir build; cd build
#export NVCC_WRAPPER_EXTRA_ARGS="-lineinfo -Xptxas=-v"
#export NVCC_WRAPPER_EXTRA_ARGS="--maxrregcount=64"
#echo $NVCC_WRAPPER_EXTRA_ARGS
#cmake .. -DCMAKE_BUILD_TYPE=Release          -DKokkos_ENABLE_CUDA=ON   -DKokkos_ARCH_VOLTA70=ON    -DKokkos_ENABLE_SERIAL=ON
cmake .. -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_VOLTA70=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_FORCE_DEVICE=ON  #-DKokkos_ENABLE_HWLOC=OFF
make #VERBOSE=1