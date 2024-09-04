#include <iostream>
#include "Lattice.h"
#include "Model.h"
#include "MonteCarlo.h"

#include <chrono>
using namespace std::chrono;

#include "ExactEnumeration.h"
using namespace std;
int main(int argc, char *argv[]) {

    //std::cout << "Hello, World!" << std::endl;
    int L = std::atoi(argv[1]);
    auto start = high_resolution_clock::now();
    std::string outFile = argv[3];
    double J = std::stod(argv[2]);
    MC_Interacting_SAW_XY mcxysaw(L, outFile);
    mcxysaw.run_simulation(J);
    mcxysaw.~MC_Interacting_SAW_XY();

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<seconds>(stop - start);

    cout << "Time " << duration.count() << endl;
    return 0;
    //std::cout << "Hello, World!" << std::endl;

   /* auto start = high_resolution_clock::now();

    //vector<vector<tuple<int, int>>> conformations = get_all_conformations(10);
    //std::cout << "Number of SAWs : " << conformations.size() << std::endl;
    int L = std::atoi(argv[1]);
    MeanValues(L);

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<seconds>(stop - start);*/

    return 0;
}
