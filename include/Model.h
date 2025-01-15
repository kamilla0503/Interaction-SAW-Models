//
// Created by Kamilla Faizullina on 08.04.2024.
//

#ifndef INTERACTION_SAW_MODELS_MODEL_H
#define INTERACTION_SAW_MODELS_MODEL_H

#include "Lattice.h"
#include "observable.h"
#include<Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#ifndef OUT_Length
#define OUT_Length 4
#endif
//used to define lattice nodes without spins (SAW does not go over this node)
#ifndef NO_SAW_NODE
#define NO_SAW_NODE -1
#endif
//used to define lattice nodes without XY spins
#ifndef NO_XY_SPIN
#define NO_XY_SPIN -5
#endif

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
    //double E;

    Kokkos::View<double*, Kokkos::CudaSpace> E;
    Kokkos::View<double*, Kokkos::CudaSpace> newE;

    Kokkos::View<coord_t*, Kokkos::CudaSpace> start_conformation;
    Kokkos::View<coord_t*, Kokkos::CudaSpace> end_conformation;

    //long end_conformation;
    //long start_conformation;
    long L;
    long lattice_side_device;
    // Random number generator pool
   Kokkos::Random_XorShift64_Pool <Kokkos::Cuda> rand_pool;


    Kokkos::View<double*, Kokkos::CudaSpace> oldspin;
    Kokkos::View<coord_t*, Kokkos::CudaSpace> save_start_conformation;
    Kokkos::View<coord_t*, Kokkos::CudaSpace> save_end_conformation;

    Kokkos::View<coord_t*, Kokkos::CudaSpace> start_index_in_nodes_position;

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
    KOKKOS_FUNCTION virtual void Energy () = 0;
    //KOKKOS_FUNCTION virtual double Energy_Add_Start () = 0;
    //KOKKOS_FUNCTION virtual double Energy_Add_End () = 0;

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
    XY_SAW_LongInteraction(long length, double J);

    KOKKOS_INLINE_FUNCTION void Reconnect(short direction); //Only Geometry changes --- the same for all SAW Models


    KOKKOS_INLINE_FUNCTION void FlipMove_AddEnd (long direction, double spinValue) override;
    KOKKOS_INLINE_FUNCTION void FlipMove_AddStart(long direction, double spinValue) override;
   // KOKKOS_INLINE_FUNCTION void ClusterStep (double flipdirection);

    void SequenceOnLatticeInitialization();
    void StartConfiguration();


    void out_MC_data(std::fstream& out, long long n_steps);
    void  updateData();

//protected:
    std::valarray<bool> used_coords;

    KOKKOS_FUNCTION void Energy ();
    //KOKKOS_FUNCTION double Energy_Add_Start () ;
    //KOKKOS_FUNCTION double Energy_Add_End () ;


    mc_stats::ScalarObservable<double> energy;
    mc_stats::ScalarObservable<double> energy_2;
    mc_stats::ScalarObservable<double> energy_4;

    mc_stats::ScalarObservable<double> mags_sin;
    mc_stats::ScalarObservable<double> mags_cos;
    mc_stats::ScalarObservable<double> magnetization_2;
    mc_stats::ScalarObservable<double> magnetization_4;



};

/*
KOKKOS_INLINE_FUNCTION
void FlipMove_AddEnd_Device(FlipMoveData& data, long direction, double spinValue) {
    printf("FlipMove_AddEnd \n");
    coord_t new_point = data.map_of_contacts_int(data.ndim2 * data.end_conformation + direction);
    //coord_t new_point = flip_data.map_of_contacts_int(lattice->ndim2() * end_conformation + direction);
    //std::cout << "new_point " << new_point << std::endl;
    printf("FlipMove_AddEnd new point = %ld \n", new_point);
    double oldspin = data.sequence_on_lattice(data.start_conformation);
    //std::cout << "oldspin " << oldspin << std::endl;
    printf("FlipMove_AddEnd new point = %f \n", oldspin);
    //self-avoidance condition:
    if (data.sequence_on_lattice(new_point) != NO_XY_SPIN) return;
    //std::cout << "sequence_on_lattice(new_point)" << sequence_on_lattice(new_point) << std::endl;
    coord_t save_start_conformation;

    // delete the beginning of SAW
    save_start_conformation = data.start_conformation;
    printf("FlipMove_AddEnd save start = %ld \n", save_start_conformation);
    data.start_conformation = data.next_monomers(data.start_conformation);
    printf("FlipMove_AddEnd start = %ld \n", data.start_conformation);
    data.next_monomers(save_start_conformation) = NO_SAW_NODE;
    data.previous_monomers(data.start_conformation) = NO_SAW_NODE;
    data.sequence_on_lattice(save_start_conformation) = NO_XY_SPIN;

    //add the new monomer at the end of SAW
    data.next_monomers(data.end_conformation) = new_point;
    data.sequence_on_lattice(new_point) = spinValue; //new spin value
    data.previous_monomers(new_point) = data.end_conformation;
    data.end_conformation = new_point;
    //std::cout << "Movement of positions  " << start_conformation << " " << spinValue << std::endl;

    for (int i = 1; i < data.L; i++) {
        data.lattice_nodes_positions(i - 1) = data.lattice_nodes_positions(i);
    }
    data.lattice_nodes_positions(data.L - 1) = data.end_conformation;

    double new_E = data.E;  // Energy();

    double p1 = exp(-(data.J * (new_E - data.E)));
    double p_metropolis = Kokkos::min(1.0, p1);

    auto rand_gen = data.rand_pool->get_state();
    // Generate a random number between 0.0 and 1.0
    double q_ifaccept = rand_gen.drand(0., 1.);
    if (q_ifaccept < p_metropolis) { // accept the new state
        data.E = new_E;
        data.sequence_on_lattice(save_start_conformation) = NO_XY_SPIN;
        data.directions(save_start_conformation) = NO_SAW_NODE;
        data.directions(data.previous_monomers(data.end_conformation)) = direction;
    } else {
        //reject new state
        //delete end
        coord_t del = data.end_conformation;
        data.end_conformation = data.previous_monomers(data.end_conformation);
        data.next_monomers(data.end_conformation) = NO_SAW_NODE;
        data.previous_monomers(del) = NO_SAW_NODE;
        data.sequence_on_lattice(del) = NO_XY_SPIN;

        //add the previous beginning
        data.previous_monomers(data.start_conformation) = save_start_conformation;
        data.next_monomers(save_start_conformation) = data.start_conformation;
        data.start_conformation = save_start_conformation;
        data.sequence_on_lattice(data.start_conformation) = oldspin;

        for (int i = data.L - 1; i > 0; i--) {
            data.lattice_nodes_positions(i) = data.lattice_nodes_positions(i - 1);
        }
        data.lattice_nodes_positions(0) = data.start_conformation;
    }
    data.rand_pool->free_state(rand_gen);
}*/

#endif //INTERACTION_SAW_MODELS_MODEL_H
