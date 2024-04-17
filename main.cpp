#include <iostream>
#include "Lattice.h"
#include "Model.h"
#include "MonteCarlo.h"

#include <chrono>
using namespace std::chrono;

#include "ExactEnumeration.h"
using namespace std;
int main(int argc, char *argv[]) {

    Lattice_2D l;
    Lattice_3D p(3);
    XY_SAW_LongInteraction xysaw;
    //XY_SAW_LongInteraction xysaw2(10);
    MC_Interacting_SAW_XY mcxysaw(8);
    mcxysaw.run_simulation();
    //std::cout << "Hello, World!" << std::endl;

    auto start = high_resolution_clock::now();

    //vector<vector<tuple<int, int>>> conformations = get_all_conformations(10);
    //std::cout << "Number of SAWs : " << conformations.size() << std::endl;
    int L = std::atoi(argv[1]);
    MeanValues(L);

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<seconds>(stop - start);

// To get the value of duration use the count()
// member function on the duration object
    //cout << "Time " << duration.count() << endl;
    return 0;
}
