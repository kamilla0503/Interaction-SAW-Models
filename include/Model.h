//
// Created by Kamilla Faizullina on 08.04.2024.
//

#ifndef INTERACTION_SAW_MODELS_MODEL_H
#define INTERACTION_SAW_MODELS_MODEL_H

#include "Lattice.h"

const double PI = std::atan(1.0)*4;

class Model {
public:
    //Model() {};
    //Model (long length);
    inline long number_of_spins() {return L;}
    inline short ndim2() {
        short dim2 = -1;
        (lattice!= nullptr) ? dim2 = lattice->ndim2() : dim2 = -1;
        return dim2;
    }
    Lattice *lattice = nullptr;
protected:
    //Model-specific Energy function; returns double as J is expected to be double also
    virtual double Energy () = 0;

    long L; //Length of the model chain
    double E; //current value for energy; double as J
    double J = 0; //Interaction Energy
};

// Abstract Class for geometry related work
template<class SpinType>
class SAW_model : public Model {
public:
    SAW_model() {};
    SAW_model(long length);

    void Reconnect(short direction); //Only Geometry changes --- the same for all SAW Models

    virtual void FlipMove_AddEnd (short direction, SpinType spinvalue) = 0; //depends on spin variables
    virtual void FlipMove_AddStart (short direction, SpinType spinvalue) = 0; //depends on spin variables
    virtual void ClusterStep () = 0; //depends on spin variables

    void LatticeInitialization();
protected:
    std::valarray<SpinType> sequence_on_lattice;
    std::valarray<long> next_monomers;
    std::valarray<long> previous_monomers;
    long end_conformation = 0;
    long start_conformation = 0;
    std::valarray<short> directions; // n-1 edges of SAW on the lattice; //directions enumerated from o to dim2()
};

//Class for XY long-interacting Model on SAWs
class XY_SAW_LongInteraction : public  SAW_model<double> {
public:
    XY_SAW_LongInteraction() {};
    XY_SAW_LongInteraction(long length);

    void FlipMove_AddEnd (short direction, double spinValue);
    void FlipMove_AddStart(short direction, double spinValue);
    void ClusterStep ();

    void SequenceOnLatticeInitialization();
    void StartConfiguration();
protected:
    double Energy ();
};

#endif //INTERACTION_SAW_MODELS_MODEL_H
