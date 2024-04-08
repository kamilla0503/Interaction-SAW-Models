//
// Created by Kamilla Faizullina on 08.04.2024.
//

#ifndef INTERACTION_SAW_MODELS_MODEL_H
#define INTERACTION_SAW_MODELS_MODEL_H

#include "Lattice.h"

class Model {
public:
    inline long number_of_spins() {return L;}
protected:
    virtual void Energy () = 0; //Model-specific Energy function

    Lattice *lattice;
    long int L; //Length of the model chain
    long E; //current value for energy
};

class SAW_model : public Model {
public:
};

#endif //INTERACTION_SAW_MODELS_MODEL_H
