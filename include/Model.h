//
// Created by Kamilla Faizullina on 08.04.2024.
//

#ifndef INTERACTION_SAW_MODELS_MODEL_H
#define INTERACTION_SAW_MODELS_MODEL_H

#include "Lattice.h"

class Model {
public:
    //Model() {};
    //Model (long length);
    inline long number_of_spins() {return L;}
protected:
    virtual void Energy () = 0; //Model-specific Energy function

    Lattice *lattice = nullptr;
    long L; //Length of the model chain
    long E; //current value for energy
};

// Abstract Class for geometry related work
class SAW_model : public Model {
public:
    SAW_model() {};
    SAW_model(long length);

    void Reconnect(int j); //Only Geometry changes --- the same for all SAW Models
    virtual void FlipMove () = 0; //depends on spin variables
    virtual void ClusterStep () = 0; //depends on spin variables

    void LatticeInitialization();
protected:
    std::valarray<long> next_monomers;
    std::valarray<long> previous_monomers;
    long end_conformation = 0;
    long start_conformation = 0;
    std::valarray<short> directions; // n-1 edges of SAW on the lattice; //directions enumerated from o to dim2()
};

//Class for XY long-interacting Model on SAWs
class XY_SAW_LongInteraction : public  SAW_model {
public:
    XY_SAW_LongInteraction() {};
    XY_SAW_LongInteraction(long length);

    void FlipMove ();
    void ClusterStep ();

    void SequenceOnLatticeInitialization();
    void StartConfiguration();
protected:
    void Energy ();
    std::valarray<double> sequence_on_lattice;
};

#endif //INTERACTION_SAW_MODELS_MODEL_H
