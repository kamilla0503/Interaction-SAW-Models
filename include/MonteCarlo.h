//
// Created by Kamilla Faizullina on 08.04.2024.
//

#ifndef INTERACTION_SAW_MODELS_MONTECARLO_H
#define INTERACTION_SAW_MODELS_MONTECARLO_H

#include "Model.h"
#include <random>
#include <chrono>
#include <string>
#include <Kokkos_Core.hpp>

class Monte_Carlo {
public:
    //Monte_Carlo () {};
protected:
    virtual void run_simulation(float J) = 0;
};

class MC_Interacting_SAW_XY : public Monte_Carlo{
public:
    MC_Interacting_SAW_XY() {};
    MC_Interacting_SAW_XY(  int length, float J, std::string LogFile = "",
                            float Probability_Local_Update = 1.00,
                            float Probability_Reconnect = 0.05,
                            float Jmin = 0.2412f, float Jmax = 0.3062f);
    void run_simulation(float J);
protected:

    std::string LogFile;

    float p_for_local_update;
    float p_for_reconnect;

    XY_SAW_LongInteraction *model = nullptr;

};

#endif //INTERACTION_SAW_MODELS_MONTECARLO_H
