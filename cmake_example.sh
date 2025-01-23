cmake .. -DKokkos_DIR="${HOME}/kokkos/install/lib/cmake/Kokkos" -DCMAKE_BUILD_TYPE=Release
cmake .. -DCMAKE_BUILD_TYPE=Release          -DKokkos_ENABLE_CUDA=ON          -DKokkos_ARCH_AMPERE80=ON          -DKokkos_ARCH_VOLTA70=OFF          -DKokkos_ENABLE_SERIAL=OFF          -DKokkos_ENABLE_OPENMP=OFF
