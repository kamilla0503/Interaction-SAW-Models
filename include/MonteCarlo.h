//
// Created by Kamilla Faizullina on 08.04.2024.
//

#ifndef INTERACTION_SAW_MODELS_MONTECARLO_H
#define INTERACTION_SAW_MODELS_MONTECARLO_H

#include "Model.h"
#include <random>
#include <chrono>


class Monte_Carlo {
public:
    //Monte_Carlo () {};
protected:
    virtual void run_simulation(double J) = 0;
};

class MC_Interacting_SAW_XY : public Monte_Carlo{
public:
    MC_Interacting_SAW_XY() {};
    MC_Interacting_SAW_XY(  long length,
                            double Probability_Local_Update = 0.9,
                            double Probability_Reconnect = 0.1 );
    void run_simulation(double J);
protected:

    double p_for_local_update;
    double p_for_reconnect;

    XY_SAW_LongInteraction *model = nullptr;

};

#endif //INTERACTION_SAW_MODELS_MONTECARLO_H
