//
// Created by Kamilla Faizullina on 08.04.2024.
//

#include "MonteCarlo.h"
#include <iostream>
#include <fstream>

#ifndef  MC_STEPS
#define MC_STEPS 15
#endif

#ifndef URD_SEED 121
#define URD_SEED
#endif
#ifndef UID_SEED 123
#define UID_SEED
#endif

MC_Interacting_SAW_XY::MC_Interacting_SAW_XY(long length,
                                             double Probability_Local_Update,
                                             double Probability_Reconnect) {
    p_for_local_update = Probability_Local_Update;
    p_for_reconnect = p_for_local_update + Probability_Reconnect;
    model = new XY_SAW_LongInteraction(length);
}


void MC_Interacting_SAW_XY::run_simulation() {

    double mc_step_type = 0;
    short step = 0;
    double flipMoveType = 0;
    double spinvalue = 0;

    static std::uniform_real_distribution<double> distribution_urd(0.0,1.0);
    static std::mt19937 generator_urd;
    generator_urd.seed(URD_SEED);

    std::uniform_int_distribution<short> distribution_uid_steps(0, model->ndim2() - 1);
    std::mt19937 generators_steps;
    generators_steps.seed(UID_SEED);

    std::uniform_real_distribution<double> distribution_theta(0, 2.0*PI);
    std::mt19937 generators_theta;
    generators_theta.seed(URD_SEED);

    for (int i = 0; i < MC_STEPS + 20; ++i) {
        /*std::fstream myStream;
        std::string filename = "For_Debug_MapOfContacts_" + std::to_string(i)+ ".out";
        myStream.open(filename,std::fstream::out);
        myStream << "start : " << model->start_conformation << std::endl;
        myStream << "end : " << model->end_conformation << std::endl;

        for (long j =0; j< model->lattice->NumberOfNodes() ; j++){
            myStream << j << " " << model->next_monomers[j] << " " << model->sequence_on_lattice[j];
            myStream << std::endl;
        }
        myStream << std::endl;
        myStream.close();*/
        mc_step_type = distribution_urd(generator_urd) ;
        if (mc_step_type < p_for_local_update) {
            flipMoveType = distribution_urd(generator_urd) ;
            step = distribution_uid_steps(generators_steps);
            spinvalue =  distribution_theta(generators_theta);
            if (flipMoveType<0.5) {
                model->FlipMove_AddEnd(step, spinvalue);
            }
            else {
                model->FlipMove_AddStart(step, spinvalue);
            }
            //model->FlipMove();
        }
        else if (mc_step_type < p_for_reconnect) {
            //std:: cout << "Reconnect" << std::endl;
            step = distribution_uid_steps(generators_steps);
            model->Reconnect(step);
        }
        else {
            model->ClusterStep();
        }
        //if (i%10000!=0) { continue;}
    }
}