//
// Created by Kamilla Faizullina on 08.04.2024.
//

#ifndef INTERACTION_SAW_MODELS_MODEL_H
#define INTERACTION_SAW_MODELS_MODEL_H

#include "Lattice.h"
#include "observable.h"
#include<Kokkos_Core.hpp>

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
    void set_J (double J_) {J = J_;}
protected:
    //Model-specific Energy function; returns double as J is expected to be double also
    virtual double Energy () = 0;

    long L; //Length of the model chain
    double E; //current value for energy; double as J
    double J; //Interaction Energy
};

// Abstract Class for geometry related work
template<class SpinType>
class SAW_model : public Model {
public:
    SAW_model() {};
    SAW_model(long length);

    void Reconnect(short direction); //Only Geometry changes --- the same for all SAW Models

    virtual void FlipMove_AddEnd (long direction, SpinType spinvalue) = 0; //depends on spin variables
    virtual void FlipMove_AddStart (long direction, SpinType spinvalue) = 0; //depends on spin variables
    virtual void ClusterStep (double flipdirection) = 0; //depends on spin variables

    void LatticeInitialization();

    //virtual void (std::fstream& out) = 0 ;
    virtual void out_MC_data(std::fstream& out, long long n_steps) = 0 ;
    virtual void updateData() = 0;

protected:
    std::valarray<SpinType> sequence_on_lattice_h;
    Kokkos::View<SpinType*>::HostMirror h_sequence_on_lattice_h;
    std::valarray<long> next_monomers_h;
    Kokkos::View<long*> next_monomers;
    Kokkos::View<long*>  previous_monomers;
    std::valarray<long> previous_monomers_h;
    long end_conformation = 0;
    long start_conformation = 0;
    std::valarray<short> directions; // n-1 edges of SAW on the lattice; //directions enumerated from o to dim2()

    mc_stats::ScalarObservable<double> e2e_distance_2;

    long* lattice_nodes_positions_h;
    Kokkos::View<long*>::HostMirror h_lattice_nodes_positions_h;

    Kokkos::View<long*>::HostMirror h_next_monomers_h;
    Kokkos::View<long*>::HostMirror h_previous_monomers_h;

    Kokkos::View<long*> lattice_nodes_positions;
    Kokkos::View<SpinType*> sequence_on_lattice;
};

//Class for XY long-interacting Model on SAWs
class XY_SAW_LongInteraction : public  SAW_model<double> {
public:
    XY_SAW_LongInteraction() {};
    XY_SAW_LongInteraction(long length);

    void FlipMove_AddEnd (long direction, double spinValue);
    void FlipMove_AddStart(long direction, double spinValue);
    void ClusterStep (double flipdirection);

    void SequenceOnLatticeInitialization();
    void StartConfiguration();


    void out_MC_data(std::fstream& out, long long n_steps);
    void  updateData();

protected:
    std::valarray<bool> used_coords;

    double Energy ();

    mc_stats::ScalarObservable<double> energy;
    mc_stats::ScalarObservable<double> energy_2;
    mc_stats::ScalarObservable<double> energy_4;

    mc_stats::ScalarObservable<double> mags_sin;
    mc_stats::ScalarObservable<double> mags_cos;
    mc_stats::ScalarObservable<double> magnetization_2;
    mc_stats::ScalarObservable<double> magnetization_4;
};

#endif //INTERACTION_SAW_MODELS_MODEL_H
