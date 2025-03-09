//
// Created by Kamilla Faizullina on 08.04.2024.
//
#include "Model.h"
#include<iostream>
#include <random>
#include <fstream>
#include <chrono>
#include<queue>
#include <Kokkos_Random.hpp>

//used to increase length of SAWs for lattice side
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


#define URD_SEED 121

#ifndef R_POWER
#define R_POWER 3
#endif

//R_power = 3
//FIx later to avoid misunderstanding
#define exponent 1.5

template<class SpinType>
SAW_model<SpinType>::SAW_model(long length) {
    L = length;
}

//Initialize geometry
template<class SpinType>
void SAW_model<SpinType>::LatticeInitialization() {
    next_monomers_h.resize(lattice->NumberOfNodes(), NO_SAW_NODE);
    previous_monomers_h.resize(lattice->NumberOfNodes(), NO_SAW_NODE);
    directions_h.resize(lattice->NumberOfNodes(), NO_SAW_NODE); //directions enumerated from o to dim2()
    lattice_nodes_positions_h = new long[number_of_spins()]{NO_SAW_NODE};
    //lattice_nodes_positions_h.resize(number_of_spins(),NO_SAW_NODE);
}

XY_SAW_LongInteraction::XY_SAW_LongInteraction(long length, double J_) : SAW_model<double>(length) {
#ifdef REGIME_2D
    lattice = new Lattice_2D(2 * L + OUT_Length);
#else
    lattice = new Lattice_3D(0.75*L+OUT_Length);
#endif
    if (lattice != nullptr) {
        J = J_;
        LatticeInitialization();
        SequenceOnLatticeInitialization();
        StartConfiguration();
    }
    //rand_pool = Kokkos::Random_XorShift64_Pool<Kokkos::DefaultExecutionSpace>(/*seed=*/12345);
    //Kokkos::fence();
    printf("Finish all configs \n");

    //lattice_side = lattice->lattice_side;
    //rand_pool.init(12345,256);
}

void XY_SAW_LongInteraction::SequenceOnLatticeInitialization() {
    //sequence_on_lattice_h = new double [lattice->NumberOfNodes()]{NO_XY_SPIN};
    sequence_on_lattice_h.resize(lattice->NumberOfNodes(), NO_XY_SPIN);
    used_coords.resize(lattice->NumberOfNodes(), false);
}
/*
static Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> g_pool;

// A helper function to initialize it
void initializePool(unsigned int seed)
{
    g_pool = Kokkos::Random_XorShift64_Pool<Kokkos::Cuda>(seed);
}

Kokkos::Random_XorShift64_Pool<Kokkos::Cuda>& getGlobalPool()
{
    return g_pool;
}*/

void XY_SAW_LongInteraction::StartConfiguration() {

    //initializePool(12345);
    //Kokkos::View<double*> sequence_on_lattice("sequence_on_lattice", this->lattice->NumberOfNodes());
    //Kokkos::View<long*> A("A", N);
    long lattice_side_h = lattice->lattice_side ; // Assign the actual value you need here
    auto Nnodes = lattice_side_h*lattice_side_h*lattice_side_h;
    
    lattice_nodes_positions = Kokkos::View<long *, Kokkos::CudaSpace>("lattice_nodes_positions", L);
    sequence_on_lattice = Kokkos::View<double *, Kokkos::CudaSpace>("sequence_on_lattice", Nnodes);
    //Kokkos::View<double*>::HostMirror
    h_sequence_on_lattice_h = Kokkos::create_mirror_view(Kokkos::HostSpace(),sequence_on_lattice);
    h_lattice_nodes_positions_h = Kokkos::create_mirror_view(Kokkos::HostSpace(),lattice_nodes_positions);
    //Kokkos::View<double *[4][4], LayoutType, MemSpace>::HostMirror h_A = Kokkos::create_mirror_view(A);
    next_monomers = Kokkos::View<long *, Kokkos::CudaSpace>("next_monomers", Nnodes);
    previous_monomers = Kokkos::View<long *, Kokkos::CudaSpace>("previous_monomers", Nnodes);
    h_next_monomers_h = Kokkos::create_mirror_view(Kokkos::HostSpace(),next_monomers);
    h_previous_monomers_h = Kokkos::create_mirror_view(Kokkos::HostSpace(),previous_monomers);

    directions = Kokkos::View<short*, Kokkos::CudaSpace>("directions", Nnodes);
    h_directions_h = Kokkos::create_mirror_view(Kokkos::HostSpace(),directions);

#ifdef STARTDEFAULT
    start_conformation = 0;
    end_conformation = L - 1;
    lattice_nodes_positions_h[0] = start_conformation;
    lattice_nodes_positions_h[this->number_of_spins() - 1] = end_conformation;
    for (int i = 1; i < L - 1; i++) {
        previous_monomers_h[i] = i - 1;
        sequence_on_lattice_h[i] = (i%6);  //PI;
        next_monomers_h[i] = i + 1;
        lattice_nodes_positions_h[i] = i;
    }
    sequence_on_lattice_h[0] = PI;
    sequence_on_lattice_h[end_conformation] = PI; //начальная последовательность
    next_monomers_h[0] = 1;
    previous_monomers_h[L - 1] = L - 2;
    //E =  -(L-1); //Hamiltonian out of (n-1) pairs of spins
    for (int i = 0; i < L - 1; i++) {
        directions_h[i] = 0; //all directions_h are the right moves
    }
#elif STARTHALF
    coord_t middle = this->number_of_spins()/2 - 1;
    start_conformation=0;
    lattice_nodes_positions_h[0] = start_conformation;
    //First part
    int i_pos = 0;
    for (int i = 1; i < middle; i++)
    {
        previous_monomers_h[i]=lattice->map_of_contacts_int_h[lattice->ndim2()*i +1];
        sequence_on_lattice_h[i]=  (i%6);  //PI;
        next_monomers_h[i]=lattice->map_of_contacts_int_h[lattice->ndim2()*i +0];
        directions_h[i]=0;
        lattice_nodes_positions_h[i] = i;
    }
    //middle
    i_pos = middle;
    previous_monomers_h[middle]=lattice->map_of_contacts_int_h[lattice->ndim2()*middle + 1];
    sequence_on_lattice_h[middle]=PI;
    next_monomers_h[middle]=lattice->map_of_contacts_int_h[lattice->ndim2()*middle + 2];
    directions_h[middle] = 2; //Go Up
    lattice_nodes_positions_h[i_pos] = middle;
    i_pos += 1;
    middle = next_monomers_h[middle];
    previous_monomers_h[middle]=lattice->map_of_contacts_int_h[lattice->ndim2()*(middle)+ 3];
    sequence_on_lattice_h[middle]=PI;
    next_monomers_h[middle]=lattice->map_of_contacts_int_h[lattice->ndim2()*(middle) + 1];
    directions_h[middle] = 1;
    lattice_nodes_positions_h[i_pos] = middle;
    i_pos += 1;
    middle = next_monomers_h[middle];
    lattice_nodes_positions_h[i_pos] = middle;
    for (int i = this->number_of_spins()/2 + 2; i < this->number_of_spins() ; i++)
    {
        previous_monomers_h[middle ]=lattice->map_of_contacts_int_h[lattice->ndim2()*middle  +0];
        sequence_on_lattice_h[middle ]= (i%6);  //PI;
        next_monomers_h[middle ]=lattice->map_of_contacts_int_h[lattice->ndim2()*middle  +1];
        directions_h[middle] = 1;
        middle = next_monomers_h[middle];
        lattice_nodes_positions_h[i] = middle;
    }
    end_conformation=middle;
    sequence_on_lattice_h[0] = PI;
    sequence_on_lattice_h[end_conformation] = PI; //начальная последовательность
    next_monomers_h[0] = lattice->map_of_contacts_int_h[lattice->ndim2()*0 +0];;
    previous_monomers_h[end_conformation] = lattice->map_of_contacts_int_h[lattice->ndim2()*end_conformation +0]; ;
    lattice_nodes_positions_h[this->number_of_spins() - 1] = end_conformation;
    std::fstream myStream;
    std::string filename = "For_Debug_AStart.out";
    myStream.open(filename,std::fstream::out);
    for (long j =0; j < this->number_of_spins() ; j++){
       myStream << lattice_nodes_positions_h[j] << " ";
    }
    myStream << std::endl;
    for (long j =0; j< this->lattice->NumberOfNodes() ; j++){
        myStream << j << " " << next_monomers_h[j] << " " << previous_monomers_h[j] << " " <<  sequence_on_lattice_h[j];
        myStream << std::endl;
    }
    myStream << std::endl;
    myStream.close();
#endif
    for (long i = 0; i < L; ++i) {
        h_lattice_nodes_positions_h(i) = lattice_nodes_positions_h[i];  // Assuming 'raw_host_data' is a long* array
    }
    for (long i = 0; i < Nnodes; ++i) {
        h_sequence_on_lattice_h(i) = sequence_on_lattice_h[i];  // Assuming 'raw_host_data' is a long* array
    }
    for (long i = 0; i <  Nnodes; ++i) {
        h_next_monomers_h(i) = next_monomers_h[i];  // Assuming 'raw_host_data' is a long* array
    }
    for (long i = 0; i <  Nnodes; ++i) {
        h_previous_monomers_h(i) = previous_monomers_h[i];  // Assuming 'raw_host_data' is a long* array
    }
    for (long i = 0; i <  Nnodes; ++i) {
        h_directions_h(i) = directions_h[i];  // Assuming 'raw_host_data' is a long* array
    }
    Kokkos::deep_copy(lattice_nodes_positions, h_lattice_nodes_positions_h);
    Kokkos::deep_copy(sequence_on_lattice, h_sequence_on_lattice_h);
    Kokkos::deep_copy(previous_monomers, h_previous_monomers_h);
    Kokkos::deep_copy(next_monomers, h_next_monomers_h);
    Kokkos::deep_copy(directions, h_directions_h);
    lattice_side = Kokkos::View<long*, Kokkos::CudaSpace>("lattice_side", 1);
    lattice_side_host = Kokkos::create_mirror_view(Kokkos::HostSpace(),lattice_side);
    lattice_side_host(0) = lattice_side_h;
    Kokkos::deep_copy(lattice_side, lattice_side_host);
    Kokkos::fence();
      //auto lattice_nodes_positions_check = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), lattice_nodes_positions);
    std::cout << "Model creation before energy" << std::endl;

    flip_data.E = Kokkos::View<double*, Kokkos::CudaSpace>("E", 1);
    flip_data.newE = Kokkos::View<double, Kokkos::CudaSpace>("newE");
    auto E_host = Kokkos::create_mirror_view(flip_data.E);
    auto newE_host = Kokkos::create_mirror_view(flip_data.newE);
    Energy();
    Kokkos::deep_copy(E, flip_data.newE);
    Kokkos::deep_copy(newE_host, flip_data.newE); // synchronizes and copies
    E = newE_host();
    Kokkos::deep_copy(flip_data.E, E);

    //auto newE_host = Kokkos::create_mirror_view(flip_data.newE);

    //E_host(0) = E;
    //newE_host(0) =
    //Kokkos::deep_copy(flip_data.E, E_host);
    //E = Energy();

    std::cout << "Model creation after energy" << std::endl;
    //printf("Energy after all = %f \n", E);

    flip_data.spinValue = Kokkos::View<double, Kokkos::CudaSpace>("spinValue");
    flip_data.direction = Kokkos::View<long, Kokkos::CudaSpace>("direction");

    std::cout << "Model creation start two scalars" << std::endl;

    flip_data.PI = Kokkos::View<double, Kokkos::CudaSpace>("PI");
    auto PI_host = Kokkos::create_mirror_view(flip_data.PI);
    //flip_data.PI
    PI_host() = std::atan(1.0)*4;
    Kokkos::deep_copy(flip_data.PI,  PI_host);

    std::cout << "Pi preparation" << std::endl;

    flip_data.L = Kokkos::View<long, Kokkos::CudaSpace>("L");
    flip_data.lattice_side_device = Kokkos::View<long, Kokkos::CudaSpace>("lattice_side_device");
    auto L_host =  Kokkos::create_mirror_view(flip_data.L);
    auto lattice_side_host =  Kokkos::create_mirror_view(flip_data.lattice_side_device);
    L_host() = L;
    lattice_side_host() =  lattice_side_h;
    Kokkos::deep_copy(flip_data.L,  L_host);
    Kokkos::deep_copy(flip_data.lattice_side_device,  lattice_side_host);


// Copy Kokkos::View members from Lattice
    flip_data.map_of_contacts_int = lattice->map_of_contacts_int;
    flip_data.inverse_steps = lattice->inverse_steps;
    flip_data.ndim2 = lattice->ndim2();
   // printf("Finish lattice fields = %f \n", E);
// Copy Kokkos::View members from XY_SAW_LongInteraction
    flip_data.sequence_on_lattice = sequence_on_lattice;
  //  printf("Finish  seq = %f \n", E);
    flip_data.next_monomers = next_monomers;
  // printf("Finish next = %f \n", E);
    flip_data.previous_monomers = previous_monomers;
   // printf("Finish prevs = %f \n", E);
    flip_data.directions = directions;
   // printf("Finish directions = %f \n", E);
    flip_data.lattice_nodes_positions = lattice_nodes_positions;
  //  printf("Finish nodes positions = %f \n", E);

// Scalars
    flip_data.J = J;
  //  printf("Finish J = %f \n", flip_data.J);
    //flip_data.E = E;
    //printf("Finish E= %f \n", E);

    //flip_data.start_conformation = start_conformation;
    //flip_data.end_conformation = end_conformation;
    //flip_data.L = L;
    //flip_data.lattice_side_device = lattice->lattice_side;
// Random pool
    //rand_pool = Kokkos::Random_XorShift64_Pool<Kokkos::Cuda>(); //(17 /* seed or execution space */);
    rand_pool = Kokkos::Random_XorShift64_Pool<Kokkos::Cuda>();
    rand_pool.init(12345,256);
    flip_data.rand_pool = rand_pool;

    //std::cout << "Start Configuration Pool states = " << rand_pool.get_num_states() << std::endl;

    //rand_pool_host = Kokkos::Random_XorShift64_Pool<Kokkos::Cuda>(17 /* seed or execution space */);
    // Allocate views with size 1
    flip_data.start_conformation = Kokkos::View<coord_t*, Kokkos::CudaSpace>("start_conformation", 1);
    flip_data.end_conformation = Kokkos::View<coord_t*, Kokkos::CudaSpace>("end_conformation", 1);

// Create host mirrors
    auto start_conformation_host = Kokkos::create_mirror_view(flip_data.start_conformation);
    auto end_conformation_host = Kokkos::create_mirror_view(flip_data.end_conformation);

// Initialize values on the host
    start_conformation_host(0) = start_conformation; // Your initial value
    end_conformation_host(0) = end_conformation;     // Your initial value

// Copy values to device
    Kokkos::deep_copy(flip_data.start_conformation, start_conformation_host);
    Kokkos::deep_copy(flip_data.end_conformation, end_conformation_host);


    flip_data.save_start_conformation = Kokkos::View<coord_t*, Kokkos::CudaSpace>("save_start_conformation", 1);
    flip_data.oldspin = Kokkos::View<double*, Kokkos::CudaSpace>("oldspin", 1);
    auto save_start_conformation_host = Kokkos::create_mirror_view(flip_data.save_start_conformation);
    auto oldspin_host = Kokkos::create_mirror_view(flip_data.oldspin);
    save_start_conformation_host(0) = -1; // Your initial value
    oldspin_host(0) = -1000;     // Your initial value



    flip_data.save_end_conformation = Kokkos::View<coord_t*, Kokkos::CudaSpace>("save_end_conformation", 1);
    auto save_end_conformation_host = Kokkos::create_mirror_view(flip_data.save_end_conformation);
    save_end_conformation_host(0) = -1; // Your initial value

    flip_data.start_index_in_nodes_position =  Kokkos::View<coord_t*, Kokkos::CudaSpace>("start_index_in_nodes_position", 1);
    auto start_index_in_nodes_position_host = Kokkos::create_mirror_view(flip_data.start_index_in_nodes_position);
    start_index_in_nodes_position_host(0) = 0;
    Kokkos::deep_copy(flip_data.start_index_in_nodes_position, start_index_in_nodes_position_host);


    // Try add i and j pairs
    long pairs_number = L*(L-1)/2;

    flip_data.i_index = Kokkos::View<long*, Kokkos::CudaSpace>("i_index", pairs_number);
    flip_data.j_index = Kokkos::View<long*, Kokkos::CudaSpace>("j_index", pairs_number);
    auto i_index_host = Kokkos::create_mirror_view(flip_data.i_index);
    auto j_index_host = Kokkos::create_mirror_view(flip_data.j_index);


    int i_pair = 0 ;
    for (int i =0; i < L; i++) {
        for (int j = i + 1; j < L; j++) {
            i_index_host(i_pair) = i;
            j_index_host(i_pair) = j;
            i_pair += 1;
        }
    }
    Kokkos::deep_copy(flip_data.i_index, i_index_host);
    Kokkos::deep_copy(flip_data.j_index, j_index_host);

    flip_data.N_pairs = Kokkos::View<long, Kokkos::CudaSpace>("N_pairs");

    auto N_pairs_host = Kokkos::create_mirror_view(flip_data.N_pairs);
    N_pairs_host() = pairs_number;
    Kokkos::deep_copy(flip_data.N_pairs,  N_pairs_host);


    flip_data.accept_move = Kokkos::View<bool, Kokkos::CudaSpace>("accept_move"); 
    auto accept_move_host = Kokkos::create_mirror_view(flip_data.accept_move);
    accept_move_host() = pairs_number;
    Kokkos::deep_copy(flip_data.accept_move,   accept_move_host);

    flip_data.flipMoveType = Kokkos::View<double, Kokkos::CudaSpace>("flipMoveType");

}

KOKKOS_INLINE_FUNCTION
double radius(const coord_t& start, const coord_t& end, long lattice_side) {
    long start_x = start % lattice_side;
    long start_y = (start % (lattice_side * lattice_side)) /lattice_side;
    long start_z = start / (lattice_side * lattice_side);
    long end_x = end % lattice_side;
    long end_y = (end % (lattice_side * lattice_side)) /lattice_side;
    long end_z = end / (lattice_side * lattice_side);

    //torus distance;
    double xdiff = abs(end_x - start_x);
    if (xdiff > (lattice_side/2))
        xdiff = lattice_side - xdiff;

    double ydiff = abs(end_y - start_y);
    if (ydiff > (lattice_side / 2))
        ydiff = lattice_side - ydiff;

    double zdiff = abs(end_z - start_z);
    if (zdiff > (lattice_side / 2))
        zdiff = lattice_side - zdiff;

    double r = xdiff *xdiff  + ydiff*ydiff + zdiff*zdiff;

    return r;
}

//this verison works and was used !
/*
KOKKOS_FUNCTION
double XY_SAW_LongInteraction::Energy() {
    double H = 0.0;  // Total energy
    const long local_L = L;
    const double lattice_side_local = lattice_side_host(0);
    auto lattice_nodes_positions_local = lattice_nodes_positions;
    auto sequence_on_lattice_local = sequence_on_lattice;

    using team_policy = Kokkos::TeamPolicy<Kokkos::Cuda>;
    using member_type = team_policy::member_type;

    // Determine the team size (you can experiment with different values)
    const int team_size = 32;  // or Kokkos::AUTO

    // Launch the parallel_reduce with team policy
    Kokkos::parallel_reduce(
            team_policy(local_L, team_size),
            KOKKOS_LAMBDA(const member_type& team_member, double& H_total) {
        const long i = team_member.league_rank();  // Get the 'i' index

        double energy_i = 0.0;

        // Parallelize the inner loop over 'j' within the team
        Kokkos::parallel_reduce(
                Kokkos::TeamThreadRange(team_member, i + 1, local_L),
                [=](const long j, double& inner_energy) {
                    double r = radius(
                            lattice_nodes_positions_local(i),
                            lattice_nodes_positions_local(j),
                            lattice_side_local
                    );
                    r = Kokkos::exp(exponent * Kokkos::log(r)); //Kokkos::pow(r, R_POWER / 2.0);

                    inner_energy += Kokkos::cos(
                            sequence_on_lattice_local(lattice_nodes_positions_local(i)) -
                            sequence_on_lattice_local(lattice_nodes_positions_local(j))
                    ) / r;
                },
                energy_i
        );
        // Each team contributes to the total energy
        Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
            H_total += energy_i;
        });
    },H);
    return -H;  // Return negative of the total energy
}*/


KOKKOS_FUNCTION
void XY_SAW_LongInteraction::Energy() {
    double H = 0.0;  // Total energy
    const long local_L = L;
    const double lattice_side_local = lattice_side_host(0);
    auto lattice_nodes_positions_local = lattice_nodes_positions;
    auto sequence_on_lattice_local = sequence_on_lattice;
    auto flip_data_local = flip_data;
    using team_policy = Kokkos::TeamPolicy<Kokkos::Cuda>;
    using member_type = team_policy::member_type;

    // Determine the team size (you can experiment with different values)
    const int team_size = 128;  // or Kokkos::AUTO

    // Launch the parallel_reduce with team policy
    Kokkos::parallel_reduce(
            team_policy(local_L, team_size),
            KOKKOS_LAMBDA(const member_type& team_member, double& H_total) {
        const long i = team_member.league_rank();  // Get the 'i' index
        double energy_i = 0.0;
        const auto pos_i = lattice_nodes_positions_local(i);
        const double theta_i = sequence_on_lattice_local(pos_i);
        // Parallelize the inner loop over 'j' within the team
        Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange(team_member, i + 1, local_L),
                [=](const long j, double& inner_energy) {
                    const auto pos_j = lattice_nodes_positions_local(j);
                    const double theta_j = sequence_on_lattice_local(pos_j);
                    double r_val = radius(pos_i, pos_j, lattice_side_local);
                    // is it faster? is it correct?
                    r_val = Kokkos::sqrt(r_val)*Kokkos::sqrt(r_val)*Kokkos::sqrt(r_val); //Kokkos::pow(r_val, exponent); // replace exp(log()) chain with pow()
                    inner_energy += Kokkos::cos(theta_i - theta_j) / r_val;

                },
                energy_i
        );
        // Each team contributes to the total energy
        Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
            H_total -= energy_i;
        });


        /*
        H_total -= energy_i;
        // (Optionally, if you want to print only once per team you might do:)
        if (team_member.team_rank() == 0) {
            printf(" i = %ld    e_i = %f \n", i, energy_i);
        }*/

    }, flip_data_local.newE );
}

std::uniform_real_distribution<double> distribution_urd(0.0, 1.0);
#ifdef SEED
std::mt19937 generator(URD_SEED + 1);
#else
std::mt19937 generator(std::chrono::steady_clock::now().time_since_epoch().count());
#endif


//My favourite and only working version for hierarchicalEnergy
//Love it
KOKKOS_INLINE_FUNCTION
void hierarchicalEnergy1(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                        const FlipMoveData &flip_data)
{
    Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team_member,  flip_data.L() ),
            [&](const long i, double &H_total) {
                double energy_i = 0.0;
                coord_t pos_i   = flip_data.lattice_nodes_positions(i);
                double theta_i  = flip_data.sequence_on_lattice(pos_i);
                const long num_j = flip_data.L() - (i + 1);
                Kokkos::parallel_reduce(
                        Kokkos::ThreadVectorRange(team_member, num_j),
                        [=](const long jj, double &innerSum) {
                            const long j = i + jj + 1;
                            coord_t pos_j  = flip_data.lattice_nodes_positions(j);
                            double theta_j = flip_data.sequence_on_lattice(pos_j);
                            double r_val   = radius(pos_i, pos_j, flip_data.lattice_side_device());
                            // Compute r_val^1.5 as before.
                            r_val = Kokkos::sqrt(r_val) * Kokkos::sqrt(r_val) * Kokkos::sqrt(r_val);
                            innerSum += Kokkos::cos(theta_i - theta_j) / r_val;
                        },
                        energy_i
                );
                    H_total -= energy_i;
            },
            flip_data.newE()
    );
}



//Try new loop
KOKKOS_INLINE_FUNCTION
void hierarchicalEnergy(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                        const FlipMoveData &flip_data)
{
    Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team_member,  flip_data.N_pairs() ),
            [&](const long ind, double &H_total) {
                    if (flip_data.accept_move()) {
                    long i = flip_data.i_index(ind);
                    long j = flip_data.j_index(ind);
                    coord_t pos_i   = flip_data.lattice_nodes_positions(i);
                    double theta_i  = flip_data.sequence_on_lattice(pos_i);
                    coord_t pos_j  = flip_data.lattice_nodes_positions(j);
                    double theta_j = flip_data.sequence_on_lattice(pos_j);
                    double r_val   = radius(pos_i, pos_j, flip_data.lattice_side_device());
                    r_val = Kokkos::sqrt(r_val) * Kokkos::sqrt(r_val) * Kokkos::sqrt(r_val);
                    H_total -= Kokkos::cos(theta_i - theta_j) / r_val;
                }
            },
            flip_data.newE()
    );
}


KOKKOS_INLINE_FUNCTION
void hierarchicalFlipMoveAddEnd(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                                const  FlipMoveData &flip_data_local,
                                Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> pool)
{
    // We'll store a bool `accept_move`. If it's false, we skip
    //bool accept_move = true;
    // single => only 1 thread in this team does the update
    Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
        // Example random usage

        auto rand_gen =  pool.get_state(); //flip_data_local.rand_pool.get_state();
        long dir = rand_gen.urand64() % 6;
        pool.free_state(rand_gen);
        flip_data_local.direction() = dir;
       // printf("hierarchicalFlipMoveAddEnd dir  = %ld; end = %ld;   \n ",   flip_data_local.direction(),flip_data_local.end_conformation(0));

        coord_t new_point = flip_data_local.map_of_contacts_int(flip_data_local.ndim2 * flip_data_local.end_conformation(0) + dir);

        // Check self-avoid
        if (flip_data_local.sequence_on_lattice(new_point) != NO_XY_SPIN) {
            //accept_move = false;
            flip_data_local.accept_move() = 0;
            return;  // skip the rest
        }
        flip_data_local.accept_move() = 1;
        auto rand_gen1 = pool.get_state();
        flip_data_local.spinValue() = rand_gen1.drand(0, 2.0*flip_data_local.PI() );
        pool.free_state(rand_gen1);
        flip_data_local.oldspin(0) = flip_data_local.sequence_on_lattice(flip_data_local.start_conformation(0));

        // delete the beginning of SAW
        flip_data_local.save_start_conformation(0) = flip_data_local.start_conformation(0);
        flip_data_local.start_conformation(0) = flip_data_local.next_monomers(flip_data_local.start_conformation(0));
        flip_data_local.next_monomers(flip_data_local.save_start_conformation(0)) = NO_SAW_NODE;
        flip_data_local.previous_monomers(flip_data_local.start_conformation(0)) = NO_SAW_NODE;
        flip_data_local.sequence_on_lattice(flip_data_local.save_start_conformation(0)) = NO_XY_SPIN;

        //add the new monomer at the end of SAW
        flip_data_local.next_monomers(flip_data_local.end_conformation(0)) = new_point;
        flip_data_local.sequence_on_lattice(new_point) = flip_data_local.spinValue(); //new spin value
        flip_data_local.previous_monomers(new_point) = flip_data_local.end_conformation(0);
        flip_data_local.end_conformation(0) = new_point;

        long position_new = flip_data_local.start_index_in_nodes_position(0) ;

        flip_data_local.lattice_nodes_positions(position_new) = flip_data_local.end_conformation(0);

    });

    // barrier if you need all threads to see the updated structure
 //   team_member.team_barrier();

    //return accept_move;
}


KOKKOS_INLINE_FUNCTION
void hierarchicalFlipMoveAddStart(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                                  const  FlipMoveData &flip_data_local,
                                  Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> pool)
{
    Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
        // Example random usage
        auto rand_gen =  pool.get_state(); //flip_data_local.rand_pool.get_state();
        flip_data_local.direction()  = rand_gen.urand64() % 6;
        pool.free_state(rand_gen);

        coord_t new_point = flip_data_local.map_of_contacts_int(flip_data_local.ndim2 * flip_data_local.start_conformation(0) + flip_data_local.direction() );
        // printf("FlipMove_AddStart new_point = %ld; start = %ld \n ",  new_point,flip_data_local.start_conformation(0));
        flip_data_local.oldspin(0) = flip_data_local.sequence_on_lattice(flip_data_local.end_conformation(0));

        if (flip_data_local.sequence_on_lattice(new_point) != NO_XY_SPIN)  {
            flip_data_local.accept_move() = 0; // Set the flag to indicate rejection
            return;
        }
        flip_data_local.accept_move() = 1;
        //coord_t flip_data_local.save_end_conformation(0);
        auto rand_gen1 = pool.get_state();
        flip_data_local.spinValue() = rand_gen1.drand(0, 2.0*flip_data_local.PI());
        pool.free_state(rand_gen1);
        //delete end
        flip_data_local.save_end_conformation(0) = flip_data_local.end_conformation(0);
        flip_data_local.end_conformation(0) = flip_data_local.previous_monomers( flip_data_local.end_conformation(0));
        flip_data_local.previous_monomers(flip_data_local.save_end_conformation(0)) = NO_SAW_NODE;
        flip_data_local.next_monomers( flip_data_local.end_conformation(0)) = NO_SAW_NODE;
        flip_data_local.sequence_on_lattice(flip_data_local.save_end_conformation(0)) = NO_XY_SPIN;

        //add the new beginning
        flip_data_local.previous_monomers(flip_data_local.start_conformation(0)) = new_point;
        flip_data_local.sequence_on_lattice(new_point) = flip_data_local.spinValue(); //выбор спина
        flip_data_local.next_monomers(new_point) = flip_data_local.start_conformation(0);
        flip_data_local.start_conformation(0) = new_point;

        long position_new = (flip_data_local.start_index_in_nodes_position(0) + flip_data_local.L() - 1) % flip_data_local.L() ;
        //if (position_new == -1 ) position_new = flip_data_local.L - 1;
        flip_data_local.lattice_nodes_positions(position_new) = flip_data_local.start_conformation(0);

    });

    // barrier if you need all threads to see the updated structure
    //   team_member.team_barrier();

   // return accept_move;
}


KOKKOS_INLINE_FUNCTION
void hierarchicalOneKernel_AddStart_FirstPart(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                                              const  FlipMoveData &flip_data_local,
                                              Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> pool)
{
    // 1) Attempt move
    //hierarchicalFlipMoveAddStart(team_member, flip_data_local, pool);
    //printf("hierarchicalOneKernel dir  = %ld; end = %ld;   \n ",   flip_data_local.direction(),flip_data_local.end_conformation(0));

    // team_member.team_barrier(); Do I need it?

    //hierarchicalEnergy(team_member, flip_data_local);
    if (!flip_data_local.accept_move() ) {
        return;
    }
   // hierarchicalEnergy(team_member, flip_data_local);

    Kokkos::single(Kokkos::PerTeam(team_member), [&]() {

       // printf("hierarchicalOneKernel Add Start Energy %f \n", flip_data_local.newE());
        double p1 = exp(-(flip_data_local.J * (flip_data_local.newE() - flip_data_local.E(0))));
        double p_metropolis = Kokkos::min(1.0, p1);
        auto rand_gen = pool.get_state();
        double q_ifaccept = rand_gen.drand(0., 1.);
        pool.free_state(rand_gen);
        if (q_ifaccept < p_metropolis) {
            flip_data_local.E(0) = flip_data_local.newE();
            flip_data_local.sequence_on_lattice(flip_data_local.save_end_conformation(0)) = NO_XY_SPIN;
            flip_data_local.directions(flip_data_local.end_conformation(0)) = NO_SAW_NODE;
            flip_data_local.directions(flip_data_local.start_conformation(0)) = flip_data_local.inverse_steps(flip_data_local.direction());
            // new start is the new added value
            long position_new = (flip_data_local.start_index_in_nodes_position(0) + flip_data_local.L () - 1) % flip_data_local.L() ;
            flip_data_local.start_index_in_nodes_position(0) = position_new;

        }
        else {
            //reject the new state
            //delete starte
            coord_t del = flip_data_local.start_conformation(0);
            flip_data_local.start_conformation(0) = flip_data_local.next_monomers(flip_data_local.start_conformation(0));
            flip_data_local.previous_monomers(flip_data_local.start_conformation(0)) = NO_SAW_NODE;
            flip_data_local.next_monomers(del) = NO_SAW_NODE;
            flip_data_local.sequence_on_lattice(del) = NO_XY_SPIN;

            //readd the end of the saw
            flip_data_local.next_monomers(flip_data_local.end_conformation(0)) = flip_data_local.save_end_conformation(0);
            flip_data_local.previous_monomers(flip_data_local.save_end_conformation(0)) = flip_data_local.end_conformation(0);
            flip_data_local.end_conformation(0) = flip_data_local.save_end_conformation(0);
            flip_data_local.sequence_on_lattice(flip_data_local.end_conformation(0)) = flip_data_local.oldspin(0);

            long position_new = (flip_data_local.start_index_in_nodes_position(0) + flip_data_local.L() - 1) % flip_data_local.L() ;

            flip_data_local.lattice_nodes_positions(position_new) = flip_data_local.end_conformation(0);
        }
        flip_data_local.newE() = 0;


    });
}

KOKKOS_INLINE_FUNCTION
void hierarchicalOneKernel_AddEnd_FirstPart(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                           const  FlipMoveData &flip_data_local,
                           Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> pool)
{
    // 1) Attempt move
    //hierarchicalFlipMoveAddEnd(team_member, flip_data_local, pool);
    //printf("hierarchicalOneKernel dir  = %ld; end = %ld;   \n ",   flip_data_local.direction(),flip_data_local.end_conformation(0));
   // team_member.team_barrier(); Do I need it?

    //hierarchicalEnergy(team_member, flip_data_local);

    if (!flip_data_local.accept_move()) {
        return;
    }
    //hierarchicalEnergy(team_member, flip_data_local);
   // team_member.team_barrier();

    Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
       // printf("hierarchicalOneKernel Energy %f \n", flip_data_local.newE());
        //printf("single dir  = %ld; end = %ld;   \n ",   flip_data_local.direction(),flip_data_local.end_conformation(0));
      //  printf("hierarchicalOneKernel Energy %f \n", flip_data_local.newE());
        double p1 = exp( -(flip_data_local.J * (flip_data_local.newE() - flip_data_local.E(0))) );
        double p_metropolis = (p1 < 1.0) ? p1 : 1.0;

        auto rand_gen = pool.get_state();
        double q_ifaccept = rand_gen.drand(0., 1.);
        pool.free_state(rand_gen);

        if (q_ifaccept < p_metropolis) {
            flip_data_local.E(0) = flip_data_local.newE();
            flip_data_local.sequence_on_lattice(flip_data_local.save_start_conformation(0)) = NO_XY_SPIN;
            flip_data_local.directions(flip_data_local.save_start_conformation(0)) = NO_SAW_NODE;
            flip_data_local.directions(flip_data_local.previous_monomers(flip_data_local.end_conformation(0))) = flip_data_local.direction();
            flip_data_local.start_index_in_nodes_position(0) = (flip_data_local.start_index_in_nodes_position(0) + 1) % flip_data_local.L();

        } else {
            // reject => revert
            // e.g. remove newly added monomer, restore old
            // ...
            coord_t del = flip_data_local.end_conformation(0);
            flip_data_local.end_conformation(0) = flip_data_local.previous_monomers(flip_data_local.end_conformation(0));
            flip_data_local.next_monomers(flip_data_local.end_conformation(0)) = NO_SAW_NODE;
            flip_data_local.previous_monomers(del) = NO_SAW_NODE;
            flip_data_local.sequence_on_lattice(del) = NO_XY_SPIN;

            //add the previous beginning
            flip_data_local.previous_monomers(flip_data_local.start_conformation(0)) = flip_data_local.save_start_conformation(0);
            flip_data_local.next_monomers(flip_data_local.save_start_conformation(0)) = flip_data_local.start_conformation(0);
            flip_data_local.start_conformation(0) = flip_data_local.save_start_conformation(0);
            flip_data_local.sequence_on_lattice(flip_data_local.start_conformation(0)) = flip_data_local.oldspin(0);

            flip_data_local.lattice_nodes_positions(flip_data_local.start_index_in_nodes_position(0)) = flip_data_local.start_conformation(0);

        }
        flip_data_local.newE()= 0;
    });
    // optional barrier
   // team_member.team_barrier();
}


 // In your XY_SAW_LongInteraction class or wherever:
void XY_SAW_LongInteraction::runMCMCOnDevice(long long MC_STEPS=10000)
{
    // (A) Create (or re-use) a random pool only once
    static bool pool_initialized = false;
    static Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> my_pool;
    if (!pool_initialized) {
        my_pool.init(/*number of states*/ 256, /*seed*/ 12345);
        pool_initialized = true;
    }

    // (B) We'll capture a copy of flip_data (assuming it's device-accessible)
    auto flip_data_local = flip_data;
    auto local_pool = my_pool;
    // (C) We launch exactly one team, with 1023 threads, as you do now
    using team_policy = Kokkos::TeamPolicy<Kokkos::Cuda>;
    team_policy policy(1, 512, 1);
    //team_policy policy(1, 1, 512);


    // (D) Single parallel_for that spawns exactly 1 team (1 block).
    //     Inside that team, we do the entire Markov chain sequentially.
    Kokkos::parallel_for("MCMC_on_device",  policy,
      KOKKOS_LAMBDA(const team_policy::member_type &team_member)
    {

        //hierarchicalEnergy(team_member, flip_data_local); 
        // Pull one random state from the pool for this entire Markov chain:
       // Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
//        auto rand_gen = local_pool.get_state();

        // (E) The MCMC loop: sequential updates, each step depends on the last
        for (long long step = 0; step < 10000+ 20; ++step)
        {
             //double flipMoveType; 
            // Decide: AddEnd vs AddStart
            Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
                auto rand_gen = local_pool.get_state();
                flip_data_local.flipMoveType() = rand_gen.drand(0., 1.);
                local_pool.free_state(rand_gen);
                if ( flip_data_local.flipMoveType()  < 0.5) {
                    // This internally does an O(N^2) parallel_reduce for energy
                   // printf("step = %lld AddEnd \n", step);
                    hierarchicalFlipMoveAddEnd(team_member, flip_data_local, local_pool);
                    //hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool );
                } else {
                    // Same logic but for "AddStart"
                   // printf("step = %lld AddStart \n", step);
                    hierarchicalFlipMoveAddStart(team_member, flip_data_local, local_pool);
                    //hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool );
                }
                // Optional: team_member.team_barrier() if you need a sync each step
                //team_member.team_barrier();
            }); //for single block 

           // printf("step = %lld Before energy  \n", step);
            hierarchicalEnergy(team_member, flip_data_local);
           // printf("step = %lld After energy  \n", step);

            if ( flip_data_local.flipMoveType()  < 0.5) {
                // This internally does an O(N^2) parallel_reduce for energy
               // printf("step = %lld AddEnd \n", step);
                //hierarchicalFlipMoveAddEnd(team_member, flip_data_local, pool);
                hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool );
            } else {
                // Same logic but for "AddStart"
              //  printf("step = %lld AddStart \n", step);
                //hierarchicalFlipMoveAddStart(team_member, flip_data_local, pool);
                hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool );
            }

        }

        // Hand back the random state
        // }); // end for single block 
//        hierarchicalEnergy(team_member, flip_data_local);

    }); // end parallel_for

    // (F) Done! We've performed MC_STEPS sequential moves on the device,
    //     with only ONE kernel launch and no host/device sync every step.
}



//KOKKOS_INLINE_FUNCTION
void XY_SAW_LongInteraction::FlipMove_AddEnd() {
    static bool pool_initialized = false;
    static Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> my_pool;
    if (!pool_initialized) {
        my_pool.init(256, 12345);
        pool_initialized = true;
    }
    auto local_pool = my_pool;
    using team_policy = Kokkos::TeamPolicy<Kokkos::Cuda>;
    using member_type = team_policy::member_type;
    //int teamSize = 128;
   // int numTeams = (L + teamSize - 1) / teamSize;
    int vectorLength = 1;
    team_policy policy(1, 1023, 1);
    //team_policy policy(1, 32, 16); //not bad choice
    //team_policy policy(numTeams, teamSize, vectorLength);

    auto flip_data_local = flip_data;
    Kokkos::parallel_for("hierarchicalKernel", policy,
                         KOKKOS_LAMBDA(const member_type &team_member) {

        /*for (long long i = 0; i < 2; ++i) {
            //double flipMoveType = distribution_urd(generator_urd) ;
            auto rand_gen = local_pool.get_state();
            double flipMoveType = rand_gen.drand(0., 1.);
            local_pool.free_state(rand_gen);
            if (flipMoveType<0.5) {
                hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool);
                //model->FlipMove_AddEnd(step, spinvalue);
                //FlipMove_AddEnd_Device(flip_data_copy, step, spinvalue);
            }
            else {
                hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool);
                // model->FlipMove_AddStart(step, spinvalue);
            }
        } */
        hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool);
    }
    );
}


void XY_SAW_LongInteraction::FlipMove_AddStart() {
    static bool pool_initialized = false;
    static Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> my_pool;
    if (!pool_initialized) {
        my_pool.init(256, 12345);
        pool_initialized = true;
    }
    auto local_pool = my_pool;
    using team_policy = Kokkos::TeamPolicy<Kokkos::Cuda>;
    using member_type = team_policy::member_type;
    //int teamSize = 128;
    // int numTeams = (L + teamSize - 1) / teamSize;
    int vectorLength = 1;
    team_policy policy(1, 1023, 1);
    // team_policy policy(1, 32, 16); //not bad choice
    //team_policy policy(numTeams, teamSize, vectorLength);

    auto flip_data_local = flip_data;
    Kokkos::parallel_for("hierarchicalKernel", policy,
                         KOKKOS_LAMBDA(const member_type &team_member) {

                            /*
        for (long long i = 0; i < 2; ++i) {
            //double flipMoveType = distribution_urd(generator_urd) ;
            auto rand_gen = local_pool.get_state();
            double flipMoveType = rand_gen.drand(0., 1.);
            local_pool.free_state(rand_gen);


            if (flipMoveType<0.5) {
                hierarchicalOneKernel_AddEnd_FirstPart(team_member, flip_data_local, local_pool);
                //model->FlipMove_AddEnd(step, spinvalue);
                //FlipMove_AddEnd_Device(flip_data_copy, step, spinvalue);
            }
            else {
                hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool);
                // model->FlipMove_AddStart(step, spinvalue);
            }
        } */

        hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool);
        // hierarchicalOneKernel_AddStart_FirstPart(team_member, flip_data_local, local_pool);
    }
    );
}





// This is correct separated version
KOKKOS_INLINE_FUNCTION
void XY_SAW_LongInteraction::FlipMove_AddEnd1() {

    auto flip_data_local = flip_data;
    // Declare a flag variable accessible on the device
    Kokkos::View<int, Kokkos::MemoryTraits<Kokkos::Atomic>> accept_move("accept_move");
    // Initialize the flag to 1 (accept by default)
    Kokkos::deep_copy(accept_move, 1);
    Kokkos::parallel_for("FlipMove_AddEnd", 1, KOKKOS_LAMBDA(const int idx) {
        auto rand_gen = flip_data_local.rand_pool.get_state();
        flip_data_local.direction()  = rand_gen.urand64() % 6;
        flip_data_local.rand_pool.free_state(rand_gen);
        coord_t new_point = flip_data_local.map_of_contacts_int(flip_data_local.ndim2* flip_data_local.end_conformation(0) + flip_data_local.direction() );
        flip_data_local.oldspin(0) = flip_data_local.sequence_on_lattice(flip_data_local.start_conformation(0));
        if (flip_data_local.sequence_on_lattice(new_point) != NO_XY_SPIN)  {
            accept_move() = 0; // Set the flag to indicate rejection
            //flip_data_local.rand_pool.free_state(rand_gen);
            return;
        }
        auto rand_gen1 = flip_data_local.rand_pool.get_state();
        flip_data_local.spinValue() = rand_gen1.drand(0, 2.0*flip_data_local.PI() );
        flip_data_local.rand_pool.free_state(rand_gen1);
        // delete the beginning of SAW
        flip_data_local.save_start_conformation(0) = flip_data_local.start_conformation(0);
        flip_data_local.start_conformation(0) = flip_data_local.next_monomers(flip_data_local.start_conformation(0));
        flip_data_local.next_monomers(flip_data_local.save_start_conformation(0)) = NO_SAW_NODE;
        flip_data_local.previous_monomers(flip_data_local.start_conformation(0)) = NO_SAW_NODE;
        flip_data_local.sequence_on_lattice(flip_data_local.save_start_conformation(0)) = NO_XY_SPIN;

        //add the new monomer at the end of SAW
        flip_data_local.next_monomers(flip_data_local.end_conformation(0)) = new_point;
        flip_data_local.sequence_on_lattice(new_point) = flip_data_local.spinValue(); //new spin value
        flip_data_local.previous_monomers(new_point) = flip_data_local.end_conformation(0);
        flip_data_local.end_conformation(0) = new_point;

        /*
        for (int i = 1; i < flip_data_local.L ; i++) {
            flip_data_local.lattice_nodes_positions(i - 1) = flip_data_local.lattice_nodes_positions(i);
        }
        flip_data_local.lattice_nodes_positions(flip_data_local.L - 1) = flip_data_local.end_conformation(0);
*/
        //temporary replacement
        //flip_data_local.lattice_nodes_positions(0) = flip_data_local.end_conformation(0);

        //write new end to the beginning
        long position_new = flip_data_local.start_index_in_nodes_position(0) ;
        //if (position_new == -1 ) position_new = flip_data_local.L - 1;
        //flip_data_local.start_index_in_nodes_position(0) = (flip_data_local.start_index_in_nodes_position(0) + 1) % flip_data_local.L;
        flip_data_local.lattice_nodes_positions(position_new) = flip_data_local.end_conformation(0);

    });

    //Kokkos::fence();
    // Copy the flag value back to the host
    int accept_move_host = 1;
    Kokkos::deep_copy(accept_move_host, accept_move);
    // If the flag indicates rejection, exit the function
    if (accept_move_host == 0) {
        return;
    }
    //auto new_E = Energy();
    Energy();
   // Kokkos::fence();

    Kokkos::parallel_for("FlipMove_AddEnd", 1, KOKKOS_LAMBDA(const int idx) {

        double p1 = exp(-(flip_data_local.J * (flip_data_local.newE() - flip_data_local.E(0)   )));
        double p_metropolis = Kokkos::min(1.0, p1);
    
        auto rand_gen = flip_data_local.rand_pool.get_state();
        // Generate a random number between 0.0 and 1.0
        double q_ifaccept = rand_gen.drand(0., 1.);
       if (q_ifaccept < p_metropolis) { // accept the new state
           flip_data_local.E(0) = flip_data_local.newE();
           flip_data_local.sequence_on_lattice(flip_data_local.save_start_conformation(0)) = NO_XY_SPIN;
           flip_data_local.directions(flip_data_local.save_start_conformation(0)) = NO_SAW_NODE;
           flip_data_local.directions(flip_data_local.previous_monomers(flip_data_local.end_conformation(0))) = flip_data_local.direction();

           flip_data_local.start_index_in_nodes_position(0) = (flip_data_local.start_index_in_nodes_position(0) + 1) % flip_data_local.L();

           /*for (int i = 1; i < flip_data_local.L ; i++) {
               flip_data_local.lattice_nodes_positions(i - 1) = flip_data_local.lattice_nodes_positions(i);
           }
           flip_data_local.lattice_nodes_positions(flip_data_local.L - 1) = flip_data_local.end_conformation(0);
       */

       } else {
            //reject new state
            //delete end
            coord_t del = flip_data_local.end_conformation(0);
            flip_data_local.end_conformation(0) = flip_data_local.previous_monomers(flip_data_local.end_conformation(0));
            flip_data_local.next_monomers(flip_data_local.end_conformation(0)) = NO_SAW_NODE;
            flip_data_local.previous_monomers(del) = NO_SAW_NODE;
            flip_data_local.sequence_on_lattice(del) = NO_XY_SPIN;
    
            //add the previous beginning
            flip_data_local.previous_monomers(flip_data_local.start_conformation(0)) = flip_data_local.save_start_conformation(0);
            flip_data_local.next_monomers(flip_data_local.save_start_conformation(0)) = flip_data_local.start_conformation(0);
            flip_data_local.start_conformation(0) = flip_data_local.save_start_conformation(0);
            flip_data_local.sequence_on_lattice(flip_data_local.start_conformation(0)) = flip_data_local.oldspin(0);

           flip_data_local.lattice_nodes_positions(flip_data_local.start_index_in_nodes_position(0)) = flip_data_local.start_conformation(0);

            /*
            for (int i = flip_data_local.L - 1; i > 0; i--) {
                flip_data_local.lattice_nodes_positions(i) = flip_data_local.lattice_nodes_positions(i - 1);
            }
           flip_data_local.lattice_nodes_positions(0) = flip_data_local.start_conformation(0);*/
        }
        flip_data_local.newE()= 0;
       flip_data_local.rand_pool.free_state(rand_gen);
    });
    //Kokkos::fence();
}

KOKKOS_INLINE_FUNCTION
void XY_SAW_LongInteraction::FlipMove_AddStart1() {

    auto flip_data_local = flip_data;
    //double flip_data_local.flip_data_local.oldspin(0)(0);
    //coord_t flip_data_local.save_start_conformation(0);
    // Declare a flag variable accessible on the device
    Kokkos::View<int, Kokkos::MemoryTraits<Kokkos::Atomic>> accept_move("accept_move");
    // Initialize the flag to 1 (accept by default)
    Kokkos::deep_copy(accept_move, 1);

    Kokkos::parallel_for("FlipMove_AddStart", 1, KOKKOS_LAMBDA(const int idx) {
        auto rand_gen = flip_data_local.rand_pool.get_state();
        flip_data_local.direction()  = rand_gen.urand64() % 6;
       // printf("FlipMove_AddStart dir  = %ld; end = %ld;   \n ",   flip_data_local.direction(),flip_data_local.end_conformation(0));
        flip_data_local.rand_pool.free_state(rand_gen);

        coord_t new_point = flip_data_local.map_of_contacts_int(flip_data_local.ndim2 * flip_data_local.start_conformation(0) + flip_data_local.direction() );
       // printf("FlipMove_AddStart new_point = %ld; start = %ld \n ",  new_point,flip_data_local.start_conformation(0));
        flip_data_local.oldspin(0) = flip_data_local.sequence_on_lattice(flip_data_local.end_conformation(0));
    
        if (flip_data_local.sequence_on_lattice(new_point) != NO_XY_SPIN)  {
            accept_move() = 0; // Set the flag to indicate rejection
            return;
        }
        //coord_t flip_data_local.save_end_conformation(0);
        auto rand_gen1 = flip_data_local.rand_pool.get_state();
        flip_data_local.spinValue() = rand_gen1.drand(0, 2.0*flip_data_local.PI());
        flip_data_local.rand_pool.free_state(rand_gen1);
        //delete end
        flip_data_local.save_end_conformation(0) = flip_data_local.end_conformation(0);
        flip_data_local.end_conformation(0) = flip_data_local.previous_monomers( flip_data_local.end_conformation(0));
        flip_data_local.previous_monomers(flip_data_local.save_end_conformation(0)) = NO_SAW_NODE;
        flip_data_local.next_monomers( flip_data_local.end_conformation(0)) = NO_SAW_NODE;
        flip_data_local.sequence_on_lattice(flip_data_local.save_end_conformation(0)) = NO_XY_SPIN;
    
        //add the new beginning
        flip_data_local.previous_monomers(flip_data_local.start_conformation(0)) = new_point;
        flip_data_local.sequence_on_lattice(new_point) = flip_data_local.spinValue(); //выбор спина
        flip_data_local.next_monomers(new_point) = flip_data_local.start_conformation(0);
        flip_data_local.start_conformation(0) = new_point;

        long position_new = (flip_data_local.start_index_in_nodes_position(0) + flip_data_local.L() - 1) % flip_data_local.L() ;
        //if (position_new == -1 ) position_new = flip_data_local.L - 1;
        flip_data_local.lattice_nodes_positions(position_new) = flip_data_local.start_conformation(0);
    });


    //Kokkos::fence();

    // Copy the flag value back to the host
    int accept_move_host = 1;
    Kokkos::deep_copy(accept_move_host, accept_move);

    // If the flag indicates rejection, exit the function
    if (accept_move_host == 0) {
        return;
    }

    //auto  new_E = Energy();
    Energy();
  //  Kokkos::fence();
    Kokkos::parallel_for("FlipMove_AddStart", 1, KOKKOS_LAMBDA(const int idx) {
       //printf("finish Energy Add Start %f \n", flip_data_local.newE());
        double p1 = exp(-(flip_data_local.J * (flip_data_local.newE() - flip_data_local.E(0))));
        double p_metropolis = Kokkos::min(1.0, p1);
        auto rand_gen = flip_data_local.rand_pool.get_state();
        double q_ifaccept = rand_gen.drand(0., 1.);
        if (q_ifaccept < p_metropolis) {
            flip_data_local.E(0) = flip_data_local.newE();
            flip_data_local.sequence_on_lattice(flip_data_local.save_end_conformation(0)) = NO_XY_SPIN;
            flip_data_local.directions(flip_data_local.end_conformation(0)) = NO_SAW_NODE;
            flip_data_local.directions(flip_data_local.start_conformation(0)) = flip_data_local.inverse_steps(flip_data_local.direction());
            // new start is the new added value
            long position_new = (flip_data_local.start_index_in_nodes_position(0) + flip_data_local.L () - 1) % flip_data_local.L() ;
            flip_data_local.start_index_in_nodes_position(0) = position_new;

        }
        else {
            //reject the new state
            //delete starte
            coord_t del = flip_data_local.start_conformation(0);
            flip_data_local.start_conformation(0) = flip_data_local.next_monomers(flip_data_local.start_conformation(0));
            flip_data_local.previous_monomers(flip_data_local.start_conformation(0)) = NO_SAW_NODE;
            flip_data_local.next_monomers(del) = NO_SAW_NODE;
            flip_data_local.sequence_on_lattice(del) = NO_XY_SPIN;
    
            //readd the end of the saw
            flip_data_local.next_monomers(flip_data_local.end_conformation(0)) = flip_data_local.save_end_conformation(0);
            flip_data_local.previous_monomers(flip_data_local.save_end_conformation(0)) = flip_data_local.end_conformation(0);
            flip_data_local.end_conformation(0) = flip_data_local.save_end_conformation(0);
            flip_data_local.sequence_on_lattice(flip_data_local.end_conformation(0)) = flip_data_local.oldspin(0);

            long position_new = (flip_data_local.start_index_in_nodes_position(0) + flip_data_local.L() - 1) % flip_data_local.L() ;

            flip_data_local.lattice_nodes_positions(position_new) = flip_data_local.end_conformation(0);
        }
        flip_data_local.newE() = 0;
        flip_data_local.rand_pool.free_state(rand_gen);
        });

    //Kokkos::fence();
}

//KOKKOS_INLINE_FUNCTION
//template<>
KOKKOS_INLINE_FUNCTION
void XY_SAW_LongInteraction::Reconnect(short direction) {

    long c = 0; //
    coord_t step_coord = lattice->map_of_contacts_int[lattice->ndim2() * end_conformation + direction];

    // test self avoidance condition
    if (sequence_on_lattice[step_coord] == NO_XY_SPIN ||
        next_monomers[step_coord] == NO_SAW_NODE ||
        step_coord == previous_monomers[end_conformation]) {
        return;
    }

    long new_end = next_monomers[step_coord];
    next_monomers[step_coord] = end_conformation;
    directions[step_coord] = lattice->inverse_steps[direction];
    c = end_conformation;
    long int new_c;
    while (c != new_end) {
        new_c = previous_monomers[c];
        next_monomers[c] = previous_monomers[c];
        directions[c] = lattice->inverse_steps[directions[new_c]];
        c = new_c;
    }
    long int temp_prev_next = next_monomers[new_end];
    previous_monomers[end_conformation] = step_coord;
    c = end_conformation;
    while (c != new_end) {
        new_c = next_monomers[c];
        previous_monomers[new_c] = c;
        c = new_c;
    }
    end_conformation = new_end;
    previous_monomers[new_end] = temp_prev_next;
    next_monomers[new_end] = NO_SAW_NODE;
    directions[new_end] = NO_SAW_NODE;

    lattice_nodes_positions[0] = start_conformation;
    c = next_monomers[start_conformation];
    for (int i = 1; i < number_of_spins(); i++) {
        lattice_nodes_positions[i] = c;
        c = next_monomers[c];
    }

}

void XY_SAW_LongInteraction::updateData() {

    auto start_host = Kokkos::create_mirror_view(flip_data.start_conformation);
    Kokkos::deep_copy(start_host, flip_data.start_conformation);
    start_conformation = start_host(0);

    auto end_host = Kokkos::create_mirror_view(flip_data.end_conformation);
    Kokkos::deep_copy(end_host, flip_data.end_conformation);
    end_conformation = end_host(0);


    double r2 = lattice->radius(start_conformation, end_conformation);
    e2e_distance_2 << r2;

    auto E_host = Kokkos::create_mirror_view(flip_data.E);
    Kokkos::deep_copy(E_host, flip_data.E);
    E = E_host(0);

    energy << E;
    energy_2 << E * E;
    energy_4 << E * E * E * E;

    double sum_sin_1 = 0.0;
    double sum_cos_1 = 0.0;
    long int current = start_conformation;

    Kokkos::deep_copy(h_sequence_on_lattice_h, flip_data.sequence_on_lattice);
    Kokkos::deep_copy(h_lattice_nodes_positions_h, flip_data.lattice_nodes_positions);

    for (int e = 0; e < L; e++) {
        sum_sin_1 += sin(h_sequence_on_lattice_h[h_lattice_nodes_positions_h[e]]);
        sum_cos_1 += cos(h_sequence_on_lattice_h[h_lattice_nodes_positions_h[e]]);

    }

    sum_sin_1 /= L;
    sum_cos_1 /= L;

    mags_sin << sum_sin_1;
    mags_cos << sum_cos_1;

    magnetization_2 << sum_sin_1 * sum_sin_1 + sum_cos_1 * sum_cos_1;
    magnetization_4
            << (sum_sin_1 * sum_sin_1 + sum_cos_1 * sum_cos_1) * (sum_sin_1 * sum_sin_1 + sum_cos_1 * sum_cos_1);

}


void XY_SAW_LongInteraction::out_MC_data(std::fstream &out, long long n_steps) {
    out << number_of_spins() << " ";
    out << J << " ";
    out << n_steps << " ";
    out << e2e_distance_2.mean() << " " << e2e_distance_2.errorbar() << " ";
    out << energy.mean() << " " << energy.errorbar() << " ";
    out << energy_2.mean() << " " << energy_2.errorbar() << " ";
    out << energy_4.mean() << " " << energy_4.errorbar() << " ";

    out << mags_sin.mean() << " " << mags_sin.errorbar() << " ";
    out << mags_cos.mean() << " " << mags_cos.errorbar() << " ";

    out << magnetization_2.mean() << " " << magnetization_2.errorbar() << " ";
    out << magnetization_4.mean() << " " << magnetization_4.errorbar() << " ";

    out << std::endl;
}