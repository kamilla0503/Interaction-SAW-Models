//
// Created by Kamilla Faizullina on 08.04.2024.
//

#include "MonteCarlo.h"
#include <iostream>
#include <fstream>


#ifndef  MC_STEPS
#define MC_STEPS 10000000000000 //99000000 //10000000000
#endif

#define URD_SEED 121
#define UID_SEED 123

MC_Interacting_SAW_XY::MC_Interacting_SAW_XY(int length, float J, std::string LogFile_,
                                             float Probability_Local_Update,
                                             float Probability_Reconnect) {
    p_for_local_update = Probability_Local_Update;
    p_for_reconnect = p_for_local_update + Probability_Reconnect;
    model = new XY_SAW_LongInteraction(length,J);
    LogFile = LogFile_;
}


//KOKKOS_INLINE_FUNCTION
void MC_Interacting_SAW_XY::run_simulation(float J) {

    //model->set_J(J);

    std::fstream MCDataStream;
    std::string filename = LogFile + "/XY_MC" + std::to_string(model->number_of_spins()) +
           "_" + std::to_string(J) + ".out";
    MCDataStream.open(filename,std::fstream::out);

    MCDataStream << "L J MC_steps R2 R2_std E E_std E2 E2_std E4 E4_std ";
    MCDataStream << "Sin1 Sin1_std Cos1 Cos1_std Mag2 Mag2_std Mag4 Mag4_std Mag1 Mag1_std ";
    MCDataStream << "eig1 eig1_std eig2 eig2_std eig3 eig3_std ";
    MCDataStream << "R_g_2_trace R_g_2_trace_std R_g_2_direct R_g_2_direct_std asphericity asphericity_std";
    MCDataStream << std::endl;



    std::fstream defect_DataStream;
    std::string filename_defect = LogFile + "/defects_" + std::to_string(model->number_of_spins()) +
    "_" + std::to_string(J) + ".out";
    defect_DataStream.open(filename_defect,std::fstream::out);



    std::fstream angles_DataStream;
    std::string filename_angles = LogFile + "/angles_" + std::to_string(model->number_of_spins()) +
    "_" + std::to_string(J) + ".out";
    angles_DataStream.open(filename_angles,std::fstream::out);


    std::fstream dirs_DataStream;
    std::string filename_dirs = LogFile + "/dirs_" + std::to_string(model->number_of_spins()) +
    "_" + std::to_string(J) + ".out";
    dirs_DataStream.open(filename_dirs,std::fstream::out);

    // Define several window sizes (number of consecutive spins to examine)
   std::vector<int> windowSizes = {10, 15, 20, 25, 30};
   // Define several threshold values (in winding number units).
   std::vector<float> thresholds = {0.7, 0.8, 0.9, 1.0, 1.1};
   for (auto w : windowSizes) {
    for (auto t : thresholds) {
        defect_DataStream << w << "_" << t << " ";
    }
   }
   defect_DataStream << "step" << std::endl;
   
    float mc_step_type = 0;
    short step = 0;
    float flipMoveType = 0;
    float spinvalue = 0;

    std::uniform_real_distribution<float> distribution_urd(0.0,1.0);
    std::mt19937 generator_urd;

    std::uniform_int_distribution<int> distribution_uid_steps(0, model->ndim2() - 1);
    std::mt19937 generators_steps;

    std::uniform_real_distribution<float> distribution_theta(0, 2.0*PI);
    std::mt19937 generators_theta;

#ifdef SEED
    generators_theta.seed(URD_SEED);
    generators_steps.seed(UID_SEED);
    generator_urd.seed(URD_SEED);
#else
    generator_urd.seed(std::chrono::steady_clock::now().time_since_epoch().count());
    generators_steps.seed(std::chrono::steady_clock::now().time_since_epoch().count());
    generators_theta.seed(std::chrono::steady_clock::now().time_since_epoch().count());
#endif

    long long n_steps_out = 20*model->number_of_spins()*model->number_of_spins();
    long long n_steps_to_equlibrium = 400*model->number_of_spins()*model->number_of_spins();
    long long n_steps_to_update = 20*model->number_of_spins()*model->number_of_spins();

    // n_steps_to_update = 10;
    // n_steps_out = 1; 
    // n_steps_to_equlibrium = 1; 
 
    long long iters = n_steps_to_update;

    auto flip_data_copy = model->flip_data; // Capture data by value
    std:: cout << "Start MC" << std:: endl;

    for (long long i = 0; i < MC_STEPS + 20; i+=iters) {
        model->runMCMCOnDevice(n_steps_to_update);
        if (i < n_steps_to_equlibrium) continue;
        if (i%(n_steps_to_update)==0) 
        model->updateData();

        if (i%(n_steps_out)==0) {
            model->out_MC_data(MCDataStream, i);
            //model->defect(defect_DataStream, i);
            model->out_angle_data(angles_DataStream, i);
            model->out_dir_data(dirs_DataStream, i);
        }
    }

    MCDataStream.close();
}
