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

    Kokkos::initialize(Kokkos::InitializationSettings()
                               .set_disable_warnings(false)
                               .set_num_threads(1)
                               );

  //  auto start = high_resolution_clock::now();
    std::string outFile = argv[3];
    double J = std::stod(argv[2]);
    std:: cout << L << " " << J << std::endl;
    MC_Interacting_SAW_XY mcxysaw(L, outFile);
    mcxysaw.run_simulation(J);
    //mcxysaw.~MC_Interacting_SAW_XY();
    Kokkos::finalize();
    //auto stop = high_resolution_clock::now();

   // auto duration = duration_cast<seconds>(stop - start);

   // cout << "Time " << duration.count() << endl;
    return 0;

}
