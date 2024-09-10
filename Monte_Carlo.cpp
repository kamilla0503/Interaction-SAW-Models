//
// Created by Kamilla Faizullina on 08.04.2024.
//

#include "MonteCarlo.h"
#include <iostream>
#include <fstream>

#ifndef  MC_STEPS
#define MC_STEPS 10000000000
#endif

#define URD_SEED 121
#define UID_SEED 123

MC_Interacting_SAW_XY::MC_Interacting_SAW_XY(long length, std::string LogFile_,
                                             double Probability_Local_Update,
                                             double Probability_Reconnect) {
    p_for_local_update = Probability_Local_Update;
    p_for_reconnect = p_for_local_update + Probability_Reconnect;
    model = new XY_SAW_LongInteraction(length);
    LogFile = LogFile_;
}


void MC_Interacting_SAW_XY::run_simulation(double J) {

    model->set_J(J);

    std::fstream MCDataStream;
    std::string filename = LogFile + "/XY_MC" + std::to_string(model->number_of_spins()) +
           "_" + std::to_string(J) + ".out";
    MCDataStream.open(filename,std::fstream::out);

    MCDataStream << "L J MC_steps R2 R2_std E E_std E2 E2_std E4 E4_std ";
    MCDataStream << "Sin1 Sin1_std Cos1 Cos1_std Mag2 Mag2_std Mag4 Mag4_std";
    MCDataStream << std::endl;

    double mc_step_type = 0;
    short step = 0;
    double flipMoveType = 0;
    double spinvalue = 0;

    std::uniform_real_distribution<double> distribution_urd(0.0,1.0);
    std::mt19937 generator_urd;

    std::uniform_int_distribution<long> distribution_uid_steps(0, model->ndim2() - 1);
    std::mt19937 generators_steps;

    std::uniform_real_distribution<double> distribution_theta(0, 2.0*PI);
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

    long long n_steps_out = 10*model->number_of_spins()*model->number_of_spins();
    long long n_steps_to_equlibrium = 100*model->number_of_spins()*model->number_of_spins();
    long long n_steps_to_update = 10*model->number_of_spins()*model->number_of_spins();

    for (long long i = 0; i < MC_STEPS + 20; ++i) {
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
        }
        else if (mc_step_type < p_for_reconnect) {
            step = distribution_uid_steps(generators_steps);
            model->Reconnect(step);
        }
        else {
            spinvalue =  distribution_theta(generators_theta);
            model->ClusterStep(spinvalue);
        }

        if (i < n_steps_to_equlibrium) continue;

        if (i%(n_steps_to_update)==0)
        model->updateData();

        if (i%(n_steps_out)==0) {
            model->out_MC_data(MCDataStream, i);
        }
    }

    MCDataStream.close();
}
