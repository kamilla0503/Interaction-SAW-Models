#include <iostream>
#include "Lattice.h"
#include "Model.h"
#include "MonteCarlo.h"
int main() {

    Lattice_2D l;
    Lattice_3D p(3);
    XY_SAW_LongInteraction xysaw;
    //XY_SAW_LongInteraction xysaw2(10);
    MC_Interacting_SAW_XY mcxysaw(8);
    mcxysaw.run_simulation();

    std::cout << "Hello, World!" << std::endl;
    return 0;
}
