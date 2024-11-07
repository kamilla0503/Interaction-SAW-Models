//
// Created by Kamilla Faizullina on 08.04.2024.
//

#ifndef INTERACTION_SAW_MODELS_MODEL_H
#define INTERACTION_SAW_MODELS_MODEL_H

#include "Lattice.h"
#include "observable.h"
#include<Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>


const double PI = std::atan(1.0)*4;
struct FlipMoveData {
    // Device-accessible data from Lattice
    Kokkos::View<long *, Kokkos::CudaSpace> map_of_contacts_int;
    Kokkos::View<int *, Kokkos::CudaSpace> inverse_steps;
    long ndim2;

    // Device-accessible data from Model
    Kokkos::View<double *, Kokkos::CudaSpace> sequence_on_lattice;
    Kokkos::View<coord_t *, Kokkos::CudaSpace> next_monomers;
    Kokkos::View<coord_t *, Kokkos::CudaSpace> previous_monomers;
    Kokkos::View<short *, Kokkos::CudaSpace> directions;
    Kokkos::View<coord_t *, Kokkos::CudaSpace> lattice_nodes_positions;

    // Scalars
    double J;
    double E;

    // Random number generator pool
   // Kokkos::Random_XorShift64_Pool <Kokkos::Cuda>* rand_pool_ptr;
};


class Model {
public:
    //Model() {};
    //Model (long length);
    KOKKOS_INLINE_FUNCTION long number_of_spins() {return L;}
    KOKKOS_INLINE_FUNCTION short ndim2() {
        //short dim2 = -1;
        //(lattice!= nullptr) ? dim2 = lattice->ndim2() : dim2 = -1;
        return lattice->ndim2();
    }
    Lattice *lattice = nullptr;
    void set_J (double J_) {J = J_;}
//protected:
    //Model-specific Energy function; returns double as J is expected to be double also
    KOKKOS_FUNCTION virtual double Energy () = 0;

    long L; //Length of the model chain
    double E; //current value for energy; double as J
    double J; //Interaction Energy
    Kokkos::View<long*, Kokkos::CudaSpace> lattice_side; //("lattice_side", 1);
    Kokkos::View<long*, Kokkos::HostSpace>::HostMirror lattice_side_host;

    FlipMoveData flip_data;

    Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> rand_pool;

};

// Abstract Class for geometry related work
template<class SpinType>
class SAW_model : public Model {
public:
    SAW_model() {};
    SAW_model(long length);

    KOKKOS_INLINE_FUNCTION virtual void Reconnect(short direction) = 0; //Only Geometry changes --- the same for all SAW Models

    KOKKOS_INLINE_FUNCTION virtual void FlipMove_AddEnd (long direction, SpinType spinvalue) = 0; //depends on spin variables
    KOKKOS_INLINE_FUNCTION virtual void FlipMove_AddStart (long direction, SpinType spinvalue) = 0; //depends on spin variables
   // virtual void ClusterStep (double flipdirection) = 0; //depends on spin variables

    void LatticeInitialization();

    //virtual void (std::fstream& out) = 0 ;
    virtual void out_MC_data(std::fstream& out, long long n_steps) = 0 ;
    virtual void updateData() = 0;

//protected:
    std::valarray<SpinType> sequence_on_lattice_h;
    typename Kokkos::View<SpinType*, Kokkos::HostSpace>::HostMirror h_sequence_on_lattice_h;
    std::valarray<long> next_monomers_h;
    Kokkos::View<long*, Kokkos::CudaSpace> next_monomers;
    Kokkos::View<long*, Kokkos::CudaSpace>  previous_monomers;
    std::valarray<long> previous_monomers_h;
    long end_conformation = 0;
    long start_conformation = 0;
    std::valarray<short> directions_h; // n-1 edges of SAW on the lattice; //directions enumerated from o to dim2()
    Kokkos::View<short*, Kokkos::CudaSpace>  directions;

    mc_stats::ScalarObservable<double> e2e_distance_2;

    long* lattice_nodes_positions_h;
    Kokkos::View<long*, Kokkos::HostSpace>::HostMirror h_lattice_nodes_positions_h;

    Kokkos::View<long*, Kokkos::HostSpace>::HostMirror h_next_monomers_h;
    Kokkos::View<long*, Kokkos::HostSpace>::HostMirror h_previous_monomers_h;
    Kokkos::View<short*, Kokkos::HostSpace>::HostMirror h_directions_h;

    Kokkos::View<long*, Kokkos::CudaSpace> lattice_nodes_positions;
    Kokkos::View<SpinType*, Kokkos::CudaSpace> sequence_on_lattice;
};

//Class for XY long-interacting Model on SAWs
class XY_SAW_LongInteraction : public  SAW_model<double> {
public:
    XY_SAW_LongInteraction() {};
    XY_SAW_LongInteraction(long length);

    KOKKOS_INLINE_FUNCTION void Reconnect(short direction); //Only Geometry changes --- the same for all SAW Models


    KOKKOS_INLINE_FUNCTION void FlipMove_AddEnd (long direction, double spinValue);
    KOKKOS_INLINE_FUNCTION void FlipMove_AddStart(long direction, double spinValue);
   // KOKKOS_INLINE_FUNCTION void ClusterStep (double flipdirection);

    void SequenceOnLatticeInitialization();
    void StartConfiguration();


    void out_MC_data(std::fstream& out, long long n_steps);
    void  updateData();

//protected:
    std::valarray<bool> used_coords;

    KOKKOS_FUNCTION double Energy ();

    mc_stats::ScalarObservable<double> energy;
    mc_stats::ScalarObservable<double> energy_2;
    mc_stats::ScalarObservable<double> energy_4;

    mc_stats::ScalarObservable<double> mags_sin;
    mc_stats::ScalarObservable<double> mags_cos;
    mc_stats::ScalarObservable<double> magnetization_2;
    mc_stats::ScalarObservable<double> magnetization_4;



};




#endif //INTERACTION_SAW_MODELS_MODEL_H
