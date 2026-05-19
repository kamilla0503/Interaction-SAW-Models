#include <iostream>
#include "Lattice.h"
#include "Model.h"
#include "MonteCarlo.h"
#include <Kokkos_Core.hpp>

//#include <chrono>
//using namespace std::chrono;
//#ifdef KOKKOS_ENABLE_OPENMP
#define OMP_PROC_BIND spread
#include "ExactEnumeration.h"
using namespace std;
int main(int argc, char *argv[]) {

    int L = std::atoi(argv[1]);

    Kokkos::initialize(Kokkos::InitializationSettings() );

    std::string outFile = argv[3];
    float J = std::stod(argv[2]);

    float Jmin = 0.2412f;
    float Jmax = 0.3062f;
    if (argc < 3) {
 
      MC_Interacting_SAW_XY mcxysaw(L, J, outFile);
      mcxysaw.run_simulation(J);
    }
    else {
      std:: cout << L << " " << J << std::endl;
      Jmin = std::stod(argv[4]);
      Jmax = std::stod(argv[5]);
      J = -1; 
      std:: cout << L << " " << Jmin << " " << Jmax << std::endl;

      MC_Interacting_SAW_XY mcxysaw(L, J, outFile, 1.00, 0.05, Jmin, Jmax);
      mcxysaw.run_simulation(J);
    }
 

    Kokkos::finalize();

    return 0;

}
