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
SAW_model<SpinType>::SAW_model(int length) {
    L = length;
}

//Initialize geometry
template<class SpinType>
void SAW_model<SpinType>::LatticeInitialization() {
    next_monomers_h.resize(lattice->NumberOfNodes(), NO_SAW_NODE);
    previous_monomers_h.resize(lattice->NumberOfNodes(), NO_SAW_NODE);
    directions_h.resize(lattice->NumberOfNodes(), NO_SAW_NODE); //directions enumerated from o to dim2()
    lattice_nodes_positions_h = new int[number_of_spins()]{NO_SAW_NODE};
    //lattice_nodes_positions_h.resize(number_of_spins(),NO_SAW_NODE);
}

XY_SAW_LongInteraction::XY_SAW_LongInteraction(int length, float J_) : SAW_model<float>(length) {
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
    //sequence_on_lattice_h = new float [lattice->NumberOfNodes()]{NO_XY_SPIN};
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
    //Kokkos::View<float*> sequence_on_lattice("sequence_on_lattice", this->lattice->NumberOfNodes());
    //Kokkos::View<int*> A("A", N);
    int lattice_side_h = lattice->lattice_side ; // Assign the actual value you need here
    auto Nnodes = lattice_side_h*lattice_side_h*lattice_side_h;
    
    lattice_nodes_positions = Kokkos::View<int**, Kokkos::CudaSpace>("lattice_nodes_positions", N_CHAINS, L);
    sequence_on_lattice = Kokkos::View<float**, Kokkos::CudaSpace>("sequence_on_lattice", N_CHAINS, Nnodes);
    //Kokkos::View<float*>::HostMirror
    h_sequence_on_lattice_h = Kokkos::create_mirror_view(sequence_on_lattice);
    h_lattice_nodes_positions_h = Kokkos::create_mirror_view(lattice_nodes_positions);
    //Kokkos::View<float *[4][4], LayoutType, MemSpace>::HostMirror h_A = Kokkos::create_mirror_view(A);
    next_monomers = Kokkos::View<int**, Kokkos::CudaSpace>("next_monomers", N_CHAINS, Nnodes);
    previous_monomers = Kokkos::View<int**, Kokkos::CudaSpace>("previous_monomers", N_CHAINS, Nnodes);
    h_next_monomers_h = Kokkos::create_mirror_view(next_monomers);
    h_previous_monomers_h = Kokkos::create_mirror_view(previous_monomers);

    directions = Kokkos::View<short**, Kokkos::CudaSpace>("directions", N_CHAINS, Nnodes);
    h_directions_h = Kokkos::create_mirror_view(directions);

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
    for (int j =0; j < this->number_of_spins() ; j++){
       myStream << lattice_nodes_positions_h[j] << " ";
    }
    myStream << std::endl;
    for (int j =0; j< this->lattice->NumberOfNodes() ; j++){
        myStream << j << " " << next_monomers_h[j] << " " << previous_monomers_h[j] << " " <<  sequence_on_lattice_h[j];
        myStream << std::endl;
    }
    myStream << std::endl;
    myStream.close();
#endif
    for (int i = 0; i < L; ++i) {
        for (int chain = 0; chain < N_CHAINS; chain++)
        h_lattice_nodes_positions_h(chain, i) = lattice_nodes_positions_h[i];  // Assuming 'raw_host_data' is a int* array
    }
    for (int i = 0; i < Nnodes; ++i) {
        for (int chain = 0; chain < N_CHAINS; chain++)
        h_sequence_on_lattice_h(chain, i) = sequence_on_lattice_h[i];  // Assuming 'raw_host_data' is a int* array
    }
    for (int i = 0; i <  Nnodes; ++i) {
        for (int chain = 0; chain < N_CHAINS; chain++)
        h_next_monomers_h(chain, i) = next_monomers_h[i];  // Assuming 'raw_host_data' is a int* array
    }
    for (int i = 0; i <  Nnodes; ++i) {
        for (int chain = 0; chain < N_CHAINS; chain++)
        h_previous_monomers_h(chain, i) = previous_monomers_h[i];  // Assuming 'raw_host_data' is a int* array
    }
    for (int i = 0; i <  Nnodes; ++i) {
        for (int chain = 0; chain < N_CHAINS; chain++)
        h_directions_h(chain, i) = directions_h[i];  // Assuming 'raw_host_data' is a int* array
    }
    Kokkos::deep_copy(lattice_nodes_positions, h_lattice_nodes_positions_h);
    Kokkos::deep_copy(sequence_on_lattice, h_sequence_on_lattice_h);
    Kokkos::deep_copy(previous_monomers, h_previous_monomers_h);
    Kokkos::deep_copy(next_monomers, h_next_monomers_h);
    Kokkos::deep_copy(directions, h_directions_h);
    lattice_side = Kokkos::View<int*, Kokkos::CudaSpace>("lattice_side", 1);
    lattice_side_host = Kokkos::create_mirror_view(Kokkos::HostSpace(),lattice_side);
    lattice_side_host(0) = lattice_side_h;
    Kokkos::deep_copy(lattice_side, lattice_side_host);
    Kokkos::fence();
      //auto lattice_nodes_positions_check = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), lattice_nodes_positions);
    std::cout << "Model creation before energy" << std::endl;

    flip_data.E = Kokkos::View<float*, Kokkos::CudaSpace>("E", N_CHAINS);
    flip_data.newE = Kokkos::View<float*, Kokkos::CudaSpace>("newE", N_CHAINS);
    // auto E_host = Kokkos::create_mirror_view(flip_data.E);
    // auto newE_host = Kokkos::create_mirror_view(flip_data.newE);
    // //Energy();
    // Kokkos::deep_copy(E, flip_data.newE);
    // Kokkos::deep_copy(newE_host, flip_data.newE); // synchronizes and copies
    // E = newE_host();
    // Kokkos::deep_copy(flip_data.E, E);

    //auto newE_host = Kokkos::create_mirror_view(flip_data.newE);

    //E_host(0) = E;
    //newE_host(0) =
    //Kokkos::deep_copy(flip_data.E, E_host);
    //E = Energy();

    std::cout << "Model creation after energy" << std::endl;
    //printf("Energy after all = %f \n", E);

    flip_data.spinValue = Kokkos::View<float*, Kokkos::CudaSpace>("spinValue", N_CHAINS);
    flip_data.direction = Kokkos::View<int*, Kokkos::CudaSpace>("direction", N_CHAINS);

    flip_data.oldIndex = Kokkos::View<int*, Kokkos::CudaSpace>("oldIndex", N_CHAINS);
    flip_data.newIndex = Kokkos::View<int*, Kokkos::CudaSpace>("newIndex", N_CHAINS);

    std::cout << "Model creation start two scalars" << std::endl;

    flip_data.PI = Kokkos::View<float, Kokkos::CudaSpace>("PI");
    auto PI_host = Kokkos::create_mirror_view(flip_data.PI);
    //flip_data.PI
    PI_host() = std::atan(1.0)*4;
    Kokkos::deep_copy(flip_data.PI,  PI_host);

    std::cout << "Pi preparation" << std::endl;

    flip_data.L = Kokkos::View<int, Kokkos::CudaSpace>("L");
    flip_data.lattice_side_device = Kokkos::View<int, Kokkos::CudaSpace>("lattice_side_device");
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
   // flip_data.J = J;
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
    flip_data.start_conformation = Kokkos::View<coord_t*, Kokkos::CudaSpace>("start_conformation", N_CHAINS);
    flip_data.end_conformation = Kokkos::View<coord_t*, Kokkos::CudaSpace>("end_conformation", N_CHAINS);

// Create host mirrors
    auto start_conformation_host = Kokkos::create_mirror_view(flip_data.start_conformation);
    auto end_conformation_host = Kokkos::create_mirror_view(flip_data.end_conformation);

// Initialize values on the host

    for (int chain = 0; chain < N_CHAINS; chain++) {
        start_conformation_host(chain) = start_conformation; // Your initial value
        end_conformation_host(chain) = end_conformation;     // Your initial value
    }
// Copy values to device
    Kokkos::deep_copy(flip_data.start_conformation, start_conformation_host);
    Kokkos::deep_copy(flip_data.end_conformation, end_conformation_host);


    flip_data.save_start_conformation = Kokkos::View<coord_t*, Kokkos::CudaSpace>("save_start_conformation", N_CHAINS);
    flip_data.oldspin = Kokkos::View<float*, Kokkos::CudaSpace>("oldspin", N_CHAINS);
    auto save_start_conformation_host = Kokkos::create_mirror_view(flip_data.save_start_conformation);
    auto oldspin_host = Kokkos::create_mirror_view(flip_data.oldspin);
 
    for (int chain = 0; chain < N_CHAINS; chain++) {
        save_start_conformation_host(chain) = -1; // Your initial value
        oldspin_host(chain) = -1000;     // Your initial value
    }


    flip_data.save_end_conformation = Kokkos::View<coord_t*, Kokkos::CudaSpace>("save_end_conformation", N_CHAINS);
    auto save_end_conformation_host = Kokkos::create_mirror_view(flip_data.save_end_conformation);
    
    for (int chain = 0; chain < N_CHAINS; chain++) 
    save_end_conformation_host(chain) = -1; // Your initial value

    flip_data.start_index_in_nodes_position =  Kokkos::View<coord_t*, Kokkos::CudaSpace>("start_index_in_nodes_position", N_CHAINS);
    auto start_index_in_nodes_position_host = Kokkos::create_mirror_view(flip_data.start_index_in_nodes_position);
    
    for (int chain = 0; chain < N_CHAINS; chain++) 
        start_index_in_nodes_position_host(chain) = 0;
    Kokkos::deep_copy(flip_data.start_index_in_nodes_position, start_index_in_nodes_position_host);

    //Kokkos::deep_copy(flip_data.start_index_in_nodes_position, 0);

    // Try add i and j pairs
    int pairs_number = L*(L-1)/2;

    flip_data.i_index = Kokkos::View<int*, Kokkos::CudaSpace>("i_index", pairs_number);
    flip_data.j_index = Kokkos::View<int*, Kokkos::CudaSpace>("j_index", pairs_number);
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

    flip_data.N_pairs = Kokkos::View<int, Kokkos::CudaSpace>("N_pairs");

    auto N_pairs_host = Kokkos::create_mirror_view(flip_data.N_pairs);
    N_pairs_host() = pairs_number;
    Kokkos::deep_copy(flip_data.N_pairs,  N_pairs_host);


    flip_data.accept_move = Kokkos::View<int*, Kokkos::CudaSpace>("accept_move", N_CHAINS); 
    //auto accept_move_host = Kokkos::create_mirror_view(flip_data.accept_move);
    //accept_move_host() = pairs_number;
    //Kokkos::deep_copy(flip_data.accept_move,   accept_move_host);

    flip_data.flipMoveType = Kokkos::View<float*, Kokkos::CudaSpace>("flipMoveType", N_CHAINS);



   // Kokkos::deep_copy(flip_data.start_index_in_nodes_position, 0);



  //  flip_data.positions = Kokkos::View<float**>("positions", N, 3);
  //  flip_data.spins     = Kokkos::View<float**>("spins",     N, 2);
//   flip_data.localField = Kokkos::View<float**>("localField", L, 2);
//   flip_data.fieldNorm = Kokkos::View<float*>("fieldNorm",  L);



//   int relax_spins =  L ; // L / 10;
//   flip_data.spin_relax = Kokkos::View<int, Kokkos::CudaSpace>("spin_relax");
//   flip_data.chosenIndices_relax = Kokkos::View<int*, Kokkos::CudaSpace>("chosenIndices_relax", relax_spins);
//   Kokkos::deep_copy(flip_data.spin_relax, relax_spins);


  flip_data.J = Kokkos::View<float, Kokkos::CudaSpace>("J");
  Kokkos::deep_copy(  flip_data.J, J);


  flip_data.d_E_1 = Kokkos::View<float*, Kokkos::CudaSpace>("d_E_1", N_CHAINS);
  flip_data.d_E_2 = Kokkos::View<float*, Kokkos::CudaSpace>("d_E_2", N_CHAINS);

  flip_data.x_coords = Kokkos::View<int *, Kokkos::CudaSpace>("x_coords", Nnodes);
  auto x_coords_h = Kokkos::create_mirror_view(flip_data.x_coords);
  for (coord_t i = 0; i < Nnodes; i++) {
    x_coords_h[i] = i % lattice_side_h;
  }
  Kokkos::deep_copy(flip_data.x_coords, x_coords_h);



  flip_data.y_coords = Kokkos::View<int *, Kokkos::CudaSpace>("y_coords", Nnodes);
  auto y_coords_h = Kokkos::create_mirror_view(flip_data.y_coords);
  for (coord_t i = 0; i < Nnodes; i++) {
    y_coords_h[i] =  (i % (lattice_side_h * lattice_side_h)) /lattice_side_h;
  }
  Kokkos::deep_copy(flip_data.y_coords, y_coords_h);

  flip_data.z_coords = Kokkos::View<int *, Kokkos::CudaSpace>("z_coords", Nnodes);
  auto z_coords_h = Kokkos::create_mirror_view(flip_data.z_coords);
  for (coord_t i = 0; i < Nnodes; i++) {
    z_coords_h[i] = i / (lattice_side_h * lattice_side_h);
  }
  Kokkos::deep_copy(flip_data.z_coords, z_coords_h);


}

KOKKOS_INLINE_FUNCTION
float radius(const coord_t& start, const coord_t& end, int lattice_side,
    const FlipMoveData &flip_data) {
    int start_x = flip_data.x_coords[start];
    int start_y = flip_data.y_coords[start];
    int start_z = flip_data.z_coords[start];
    int end_x = flip_data.x_coords[end];
    int end_y = flip_data.y_coords[end];
    int end_z = flip_data.z_coords[end];
    //torus distance;
    float xdiff = abs(end_x - start_x);
    if (xdiff > (lattice_side/2))
        xdiff = lattice_side - xdiff;

    float ydiff = abs(end_y - start_y);
    if (ydiff > (lattice_side / 2))
        ydiff = lattice_side - ydiff;

    float zdiff = abs(end_z - start_z);
    if (zdiff > (lattice_side / 2))
        zdiff = lattice_side - zdiff;

    float r = xdiff *xdiff  + ydiff*ydiff + zdiff*zdiff;

    return r;
}

//this verison works and was used !
/*
KOKKOS_FUNCTION
float XY_SAW_LongInteraction::Energy() {
    float H = 0.0;  // Total energy
    const int local_L = L;
    const float lattice_side_local = lattice_side_host(0);
    auto lattice_nodes_positions_local = lattice_nodes_positions;
    auto sequence_on_lattice_local = sequence_on_lattice;

    using team_policy = Kokkos::TeamPolicy<Kokkos::Cuda>;
    using member_type = team_policy::member_type;

    // Determine the team size (you can experiment with different values)
    const int team_size = 32;  // or Kokkos::AUTO

    // Launch the parallel_reduce with team policy
    Kokkos::parallel_reduce(
            team_policy(local_L, team_size),
            KOKKOS_LAMBDA(const member_type& team_member, float& H_total) {
        const int i = team_member.league_rank();  // Get the 'i' index

        float energy_i = 0.0;

        // Parallelize the inner loop over 'j' within the team
        Kokkos::parallel_reduce(
                Kokkos::TeamThreadRange(team_member, i + 1, local_L),
                [=](const int j, float& inner_energy) {
                    float r = radius(
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


KOKKOS_INLINE_FUNCTION
float radius1(const coord_t& start, const coord_t& end, int lattice_side) {
    int start_x = start % lattice_side;
    int start_y = (start % (lattice_side * lattice_side)) /lattice_side;
    int start_z = start / (lattice_side * lattice_side);
    int end_x = end % lattice_side;
    int end_y = (end % (lattice_side * lattice_side)) /lattice_side;
    int end_z = end / (lattice_side * lattice_side);
    //torus distance;
    float xdiff = abs(end_x - start_x);
    if (xdiff > (lattice_side/2))
        xdiff = lattice_side - xdiff;

    float ydiff = abs(end_y - start_y);
    if (ydiff > (lattice_side / 2))
        ydiff = lattice_side - ydiff;

    float zdiff = abs(end_z - start_z);
    if (zdiff > (lattice_side / 2))
        zdiff = lattice_side - zdiff;

    float r = xdiff *xdiff  + ydiff*ydiff + zdiff*zdiff;

    return r;
}


// KOKKOS_FUNCTION
// void XY_SAW_LongInteraction::Energy() {
//     float H = 0.0;  // Total energy
//     const int local_L = L;
//     const float lattice_side_local = lattice_side_host(0);
//     auto lattice_nodes_positions_local = lattice_nodes_positions;
//     auto sequence_on_lattice_local = sequence_on_lattice;
//     auto flip_data_local = flip_data;
//     using team_policy = Kokkos::TeamPolicy<Kokkos::Cuda>;
//     using member_type = team_policy::member_type;

//     // Determine the team size (you can experiment with different values)
//     const int team_size = 128;  // or Kokkos::AUTO

//     // Launch the parallel_reduce with team policy
//     Kokkos::parallel_reduce(
//             team_policy(local_L, team_size),
//             KOKKOS_LAMBDA(const member_type& team_member, float& H_total) {
//         const int i = team_member.league_rank();  // Get the 'i' index
//         float energy_i = 0.0;
//         const auto pos_i = lattice_nodes_positions_local(i);
//         const float theta_i = sequence_on_lattice_local(pos_i);
//         // Parallelize the inner loop over 'j' within the team
//         Kokkos::parallel_reduce(
//                 Kokkos::TeamVectorRange(team_member, i + 1, local_L),
//                 [=](const int j, float& inner_energy) {
//                     const auto pos_j = lattice_nodes_positions_local(j);
//                     const float theta_j = sequence_on_lattice_local(pos_j);
//                     float r_val = radius1(pos_i, pos_j, lattice_side_local);
//                     // is it faster? is it correct?
//                     r_val = Kokkos::sqrt(r_val)*Kokkos::sqrt(r_val)*Kokkos::sqrt(r_val); //Kokkos::pow(r_val, exponent); // replace exp(log()) chain with pow()
//                     inner_energy += Kokkos::cos(theta_i - theta_j) / r_val;

//                 },
//                 energy_i
//         );
//         // Each team contributes to the total energy
//         Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
//             H_total -= energy_i;
//         });
//     }, flip_data_local.newE );
// }

std::uniform_real_distribution<float> distribution_urd(0.0, 1.0);
#ifdef SEED
std::mt19937 generator(URD_SEED + 1);
#else
std::mt19937 generator(std::chrono::steady_clock::now().time_since_epoch().count());
#endif


//My favourite and only working version for hierarchicalEnergy
//Love it
// KOKKOS_INLINE_FUNCTION
// void hierarchicalEnergy1(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
//                         const FlipMoveData &flip_data)
// {
//     Kokkos::parallel_reduce(
//             Kokkos::TeamThreadRange(team_member,  flip_data.L() ),
//             [&](const int i, float &H_total) {
//                 float energy_i = 0.0;
//                 coord_t pos_i   = flip_data.lattice_nodes_positions(i);
//                 float theta_i  = flip_data.sequence_on_lattice(pos_i);
//                 const int num_j = flip_data.L() - (i + 1);
//                 Kokkos::parallel_reduce(
//                         Kokkos::ThreadVectorRange(team_member, num_j),
//                         [=](const int jj, float &innerSum) {
//                             const int j = i + jj + 1;
//                             coord_t pos_j  = flip_data.lattice_nodes_positions(j);
//                             float theta_j = flip_data.sequence_on_lattice(pos_j);
//                             float r_val   = radius(pos_i, pos_j, flip_data.lattice_side_device(),flip_data);
//                             // Compute r_val^1.5 as before.
//                             r_val = Kokkos::sqrt(r_val) * Kokkos::sqrt(r_val) * Kokkos::sqrt(r_val);
//                             innerSum += Kokkos::cos(theta_i - theta_j) / r_val;
//                         },
//                         energy_i
//                 );
//                     H_total -= energy_i;
//             },
//             flip_data.newE()
//     );
// }

//Try new loop
KOKKOS_INLINE_FUNCTION
void hierarchicalEnergy(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                        const FlipMoveData &flip_data, int c)
{
    Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team_member,  flip_data.N_pairs() ),
            [&](const int ind, float &H_total) {
                 //   if (flip_data.accept_move()) {
                    int i = flip_data.i_index(ind);
                    int j = flip_data.j_index(ind);
                    coord_t pos_i   = flip_data.lattice_nodes_positions(c, i);
                    float theta_i  = flip_data.sequence_on_lattice(c, pos_i);
                    coord_t pos_j  = flip_data.lattice_nodes_positions(c, j);
                    float theta_j = flip_data.sequence_on_lattice(c, pos_j);
                    float r_val   = radius(pos_i, pos_j, flip_data.lattice_side_device(),flip_data);
                    r_val = Kokkos::sqrt(r_val) * Kokkos::sqrt(r_val) * Kokkos::sqrt(r_val);
                    H_total -= Kokkos::cos(theta_i - theta_j) / r_val;
                //}
            },
            flip_data.newE(c)
    );
}




//Try O(N)
//This is for old spin 
KOKKOS_INLINE_FUNCTION
void hierarchicalDeltaE_1(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                        const FlipMoveData &flip_data, int c)
{
    coord_t pos_j  = flip_data.oldIndex(c);
    float theta_j = flip_data.oldspin(c);
    coord_t pos_k  = flip_data.newIndex(c);
    float theta_k = flip_data.sequence_on_lattice(c, pos_k);
    Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team_member,  flip_data.L() ),
            [&](const int i, float &H_total) {

                coord_t pos_i = flip_data.lattice_nodes_positions(c, i);
                if ((pos_i ==  flip_data.oldIndex(c)) || (pos_i == flip_data.newIndex(c)) ) return; 

                    float theta_i  = flip_data.sequence_on_lattice(c, pos_i);
                    float r_val   = radius(pos_i, pos_j, flip_data.lattice_side_device(), flip_data);
                    r_val = Kokkos::sqrt(r_val) * r_val; //Kokkos::pow(r_val, exponent); //Kokkos::sqrt(r_val) * Kokkos::sqrt(r_val) * Kokkos::sqrt(r_val);
                    H_total += Kokkos::cos(theta_i - theta_j) / r_val;

                    r_val   = radius(pos_i, pos_k , flip_data.lattice_side_device(), flip_data);
                    r_val = Kokkos::sqrt(r_val) * r_val; //Kokkos::pow(r_val, exponent);    // Kokkos::sqrt(r_val) * Kokkos::sqrt(r_val) * Kokkos::sqrt(r_val);
                    H_total -= Kokkos::cos(theta_i - theta_k) / r_val;

            },
           flip_data.d_E_1(c)
    );
   // team_member.team_barrier();

}



KOKKOS_INLINE_FUNCTION
void hierarchicalFlipMoveAddEnd(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                               const  FlipMoveData &flip_data_local, int c, 
                               const Kokkos::Random_XorShift64_Pool<Kokkos::Cuda>& pool)
{
    // We'll store a bool `accept_move`. If it's false, we skip
    //bool accept_move = true;
    // single => only 1 thread in this team does the update
    Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
        // Example random usage

        auto rand_gen =  pool.get_state();  
        int dir = rand_gen.urand64() % 6;
        pool.free_state(rand_gen);
        flip_data_local.direction(c) = dir;

        coord_t new_point = flip_data_local.map_of_contacts_int(flip_data_local.ndim2 * flip_data_local.end_conformation(c) + dir);

        // Check self-avoid
        if (flip_data_local.sequence_on_lattice(c, new_point) != NO_XY_SPIN) {
            //accept_move = false;
            flip_data_local.accept_move(c) = 0;
            return;  // skip the rest
        }
        flip_data_local.accept_move(c) = 1;
        auto rand_gen1 = pool.get_state();
        flip_data_local.spinValue(c) = rand_gen1.drand(0, 2.0*flip_data_local.PI() );
        pool.free_state(rand_gen1);
        flip_data_local.oldspin(c) = flip_data_local.sequence_on_lattice(c, flip_data_local.start_conformation(c));
         

        // delete the beginning of SAW
        flip_data_local.save_start_conformation(c) = flip_data_local.start_conformation(c);
        flip_data_local.start_conformation(c) = flip_data_local.next_monomers(c, flip_data_local.start_conformation(c));
        flip_data_local.next_monomers(c, flip_data_local.save_start_conformation(c)) = NO_SAW_NODE;
        flip_data_local.previous_monomers(c, flip_data_local.start_conformation(c)) = NO_SAW_NODE;
        flip_data_local.sequence_on_lattice(c, flip_data_local.save_start_conformation(c)) = NO_XY_SPIN;


        flip_data_local.oldIndex(c) = flip_data_local.save_start_conformation(c);
        flip_data_local.newIndex(c) = new_point;

        //add the new monomer at the end of SAW
        flip_data_local.next_monomers(c, flip_data_local.end_conformation(c)) = new_point;
        flip_data_local.sequence_on_lattice(c, new_point) = flip_data_local.spinValue(c); //new spin value
        flip_data_local.previous_monomers(c, new_point) = flip_data_local.end_conformation(c);
        flip_data_local.end_conformation(c) = new_point;

        int position_new = flip_data_local.start_index_in_nodes_position(c) ;

        flip_data_local.lattice_nodes_positions(c, position_new) = flip_data_local.end_conformation(c);
    });
}

KOKKOS_INLINE_FUNCTION
void hierarchicalFlipMoveAddStart(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                                  const FlipMoveData &flip_data_local, int c,
                                  const Kokkos::Random_XorShift64_Pool<Kokkos::Cuda>& pool)
{
    Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
        // Example random usage
        auto rand_gen =  pool.get_state();  
        flip_data_local.direction(c)  = rand_gen.urand64() % 6;
        pool.free_state(rand_gen);

        coord_t new_point = flip_data_local.map_of_contacts_int(flip_data_local.ndim2 * flip_data_local.start_conformation(c) + flip_data_local.direction(c) );
        flip_data_local.oldspin(c) = flip_data_local.sequence_on_lattice(c, flip_data_local.end_conformation(c));

        if (flip_data_local.sequence_on_lattice(c, new_point) != NO_XY_SPIN)  {
            flip_data_local.accept_move(c) = 0; // Set the flag to indicate rejection
            return;
        }
        flip_data_local.accept_move(c) = 1;
        //coord_t flip_data_local.save_end_conformation(0);
        auto rand_gen1 = pool.get_state();
        flip_data_local.spinValue(c) = rand_gen1.drand(0, 2.0*flip_data_local.PI());
        pool.free_state(rand_gen1);
        //delete end
        flip_data_local.save_end_conformation(c) = flip_data_local.end_conformation(c);
        flip_data_local.end_conformation(c) = flip_data_local.previous_monomers(c, flip_data_local.end_conformation(c));
        flip_data_local.previous_monomers(c, flip_data_local.save_end_conformation(c)) = NO_SAW_NODE;
        flip_data_local.next_monomers(c, flip_data_local.end_conformation(c)) = NO_SAW_NODE;
        flip_data_local.sequence_on_lattice(c, flip_data_local.save_end_conformation(c)) = NO_XY_SPIN;


        flip_data_local.oldIndex(c) = flip_data_local.save_end_conformation(c);
        flip_data_local.newIndex(c) =  new_point;

        //add the new beginning
        flip_data_local.previous_monomers(c, flip_data_local.start_conformation(c)) = new_point;
        flip_data_local.sequence_on_lattice(c, new_point) = flip_data_local.spinValue(c); //выбор спина
        flip_data_local.next_monomers(c, new_point) = flip_data_local.start_conformation(c);
        flip_data_local.start_conformation(c) = new_point;

        int position_new = (flip_data_local.start_index_in_nodes_position(c) + flip_data_local.L() - 1) % flip_data_local.L() ;
        //if (position_new == -1 ) position_new = flip_data_local.L - 1;
        flip_data_local.lattice_nodes_positions(c, position_new) = flip_data_local.start_conformation(c);

    });
}


KOKKOS_INLINE_FUNCTION
void hierarchicalOneKernel_AddStart_FirstPart(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                                             const  FlipMoveData &flip_data_local, int c,
                                             const Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> & pool)
{
    if (!flip_data_local.accept_move(c) ) {
        return;
    }
    Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
        float p1 = exp(-(flip_data_local.J() * (  flip_data_local.d_E_1(c)   )));
        float p_metropolis = Kokkos::min(1.0f, p1);
        auto rand_gen = pool.get_state();
        float q_ifaccept = rand_gen.drand(0., 1.);
        pool.free_state(rand_gen);
        if (q_ifaccept < p_metropolis) {
            //flip_data_local.E(0) = flip_data_local.newE();
            flip_data_local.sequence_on_lattice(c, flip_data_local.save_end_conformation(c)) = NO_XY_SPIN;
            flip_data_local.directions(c, flip_data_local.end_conformation(c)) = NO_SAW_NODE;
            flip_data_local.directions(c, flip_data_local.start_conformation(c)) = flip_data_local.inverse_steps(flip_data_local.direction(c));
            // new start is the new added value
            int position_new = (flip_data_local.start_index_in_nodes_position(c) + flip_data_local.L () - 1) % flip_data_local.L() ;
            flip_data_local.start_index_in_nodes_position(c) = position_new;

        }
        else {
            //reject the new state
            //delete starte
            coord_t del = flip_data_local.start_conformation(c);
            flip_data_local.start_conformation(c) = flip_data_local.next_monomers(c, flip_data_local.start_conformation(c));
            flip_data_local.previous_monomers(c, flip_data_local.start_conformation(c)) = NO_SAW_NODE;
            flip_data_local.next_monomers(c, del) = NO_SAW_NODE;
            flip_data_local.sequence_on_lattice(c, del) = NO_XY_SPIN;

            //readd the end of the saw
            flip_data_local.next_monomers(c, flip_data_local.end_conformation(c)) = flip_data_local.save_end_conformation(c);
            flip_data_local.previous_monomers(c, flip_data_local.save_end_conformation(c)) = flip_data_local.end_conformation(c);
            flip_data_local.end_conformation(c) = flip_data_local.save_end_conformation(c);
            flip_data_local.sequence_on_lattice(c, flip_data_local.end_conformation(c)) = flip_data_local.oldspin(c);

            int position_new = (flip_data_local.start_index_in_nodes_position(c) + flip_data_local.L() - 1) % flip_data_local.L() ;

            flip_data_local.lattice_nodes_positions(c, position_new) = flip_data_local.end_conformation(c);
        }
        //flip_data_local.newE() = 0;


    });
}

KOKKOS_INLINE_FUNCTION
void hierarchicalOneKernel_AddEnd_FirstPart(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                        const FlipMoveData &flip_data_local, int c,
                        const Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> & pool)
{
    if (!flip_data_local.accept_move(c)) {
        return;
    }
    Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
       // float p1 = exp( -(flip_data_local.J() * (flip_data_local.newE() - flip_data_local.E(0))) );
//        float p1 = exp( -(flip_data_local.J() * ( flip_data_local.d_E_2() - flip_data_local.d_E_1()       )) );
  

float p1 = exp( -(flip_data_local.J() * (  flip_data_local.d_E_1(c)       )) );
      

        float p_metropolis = (p1 < 1.0) ? p1 : 1.0;

        auto rand_gen = pool.get_state();
        float q_ifaccept = rand_gen.drand(0., 1.);
        pool.free_state(rand_gen);

        if (q_ifaccept < p_metropolis) {
           // printf("Accept new state \n "); 
            //flip_data_local.E(0) = flip_data_local.newE();
            flip_data_local.sequence_on_lattice(c, flip_data_local.save_start_conformation(c)) = NO_XY_SPIN;
            flip_data_local.directions(c, flip_data_local.save_start_conformation(c)) = NO_SAW_NODE;
            flip_data_local.directions(c, flip_data_local.previous_monomers(c, flip_data_local.end_conformation(c))) = flip_data_local.direction(c);
            flip_data_local.start_index_in_nodes_position(c) = (flip_data_local.start_index_in_nodes_position(c) + 1) % flip_data_local.L();

        } else {
            // reject => revert
            coord_t del = flip_data_local.end_conformation(c);
            flip_data_local.end_conformation(c) = flip_data_local.previous_monomers(c, flip_data_local.end_conformation(c));
            flip_data_local.next_monomers(c, flip_data_local.end_conformation(c)) = NO_SAW_NODE;
            flip_data_local.previous_monomers(c, del) = NO_SAW_NODE;
            flip_data_local.sequence_on_lattice(c, del) = NO_XY_SPIN;

            //add the previous beginning
            flip_data_local.previous_monomers(c, flip_data_local.start_conformation(c)) = flip_data_local.save_start_conformation(c);
            flip_data_local.next_monomers(c, flip_data_local.save_start_conformation(c)) = flip_data_local.start_conformation(c);
            flip_data_local.start_conformation(c) = flip_data_local.save_start_conformation(c);
            flip_data_local.sequence_on_lattice(c, flip_data_local.start_conformation(c)) = flip_data_local.oldspin(c);

            flip_data_local.lattice_nodes_positions(c, flip_data_local.start_index_in_nodes_position(c)) = flip_data_local.start_conformation(c);

        }
        //flip_data_local.newE()= 0;
    });
}


KOKKOS_INLINE_FUNCTION
void hierarchicalOneKernel_Reconnect(const Kokkos::TeamPolicy<Kokkos::Cuda>::member_type &team_member,
                        const FlipMoveData &flip_data_local, int chain,
                        const  Kokkos::Random_XorShift64_Pool<Kokkos::Cuda>& pool) 
{
    Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
        auto rand_gen =  pool.get_state(); //flip_data_local.rand_pool.get_state();
        flip_data_local.direction(chain)  = rand_gen.urand64() % 6;
        pool.free_state(rand_gen);

        int  step_coord = flip_data_local.map_of_contacts_int(flip_data_local.ndim2 * flip_data_local.end_conformation(chain) + flip_data_local.direction(chain) );

        int c = 0; //
 
        // test self avoidance condition
        if (flip_data_local.sequence_on_lattice(chain, step_coord) == NO_XY_SPIN ||
            flip_data_local.next_monomers(chain, step_coord) == NO_SAW_NODE ||
            step_coord == flip_data_local.previous_monomers(chain, flip_data_local.end_conformation(chain))) {
            return;
        }

    int new_end = flip_data_local.next_monomers(chain, step_coord);
    flip_data_local.next_monomers(chain, step_coord) = flip_data_local.end_conformation(chain);
    //need to check inverse steps 
    flip_data_local.directions(chain, step_coord) = flip_data_local.inverse_steps(flip_data_local.direction(chain));
    c = flip_data_local.end_conformation(chain);
    int new_c;
    while (c != new_end) {
        new_c = flip_data_local.previous_monomers(chain, c);
        flip_data_local.next_monomers(chain, c) = flip_data_local.previous_monomers(chain, c);
        flip_data_local.directions(chain, c) = flip_data_local.inverse_steps(flip_data_local.directions(chain, new_c) );
        c = new_c;
    }
    int temp_prev_next = flip_data_local.next_monomers(chain, new_end);
    flip_data_local.previous_monomers(chain, flip_data_local.end_conformation(chain)) = step_coord;
    c = flip_data_local.end_conformation(chain);
    while (c != new_end) {
        new_c = flip_data_local.next_monomers(chain, c);
        flip_data_local.previous_monomers(chain, new_c) = c;
        c = new_c;
    }
    flip_data_local.end_conformation(chain) = new_end;
    flip_data_local.previous_monomers(chain, new_end) = temp_prev_next;
    flip_data_local.next_monomers(chain, new_end) = NO_SAW_NODE;
    flip_data_local.directions(chain, new_end) = NO_SAW_NODE;

    flip_data_local.lattice_nodes_positions(chain, 0) = flip_data_local.start_conformation(chain);
    c = flip_data_local.next_monomers(chain, flip_data_local.start_conformation(chain));
    for (int i = 1; i < flip_data_local.L(); i++) {
        flip_data_local.lattice_nodes_positions(chain, i) = c;
        c = flip_data_local.next_monomers(chain, c);
    }
    flip_data_local.start_index_in_nodes_position(chain) = 0;
    //Redefine positions in array now 
    });

}


 // In your XY_SAW_LongInteraction class or wherever:
void XY_SAW_LongInteraction::runMCMCOnDevice(long long MC_STEPS=10000)
{
    // (A) Create (or re-use) a random pool only once
    static bool pool_initialized = false;
    static Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> my_pool;
    if (!pool_initialized) {
        my_pool.init(/*number of states*/ 256, /*seed*/ 42);
        //my_pool.init(42);
        pool_initialized = true;
    }
    // (B) We'll capture a copy of flip_data (assuming it's device-accessible)
    auto flip_data_local = flip_data;
  // FlipMoveData flip_data_local(flip_data);
    auto pool = my_pool;
    // (C) We launch exactly one team, with 1023 threads, as you do now
    using team_policy = Kokkos::TeamPolicy<Kokkos::Cuda>;
    team_policy policy(N_CHAINS, 500, 1);
    auto n_iters = MC_STEPS;

    // (D) Single parallel_for that spawns exactly 1 team (1 block).
    //     Inside that team, we do the entire Markov chain sequentially.
    Kokkos::parallel_for("MCMC_on_device",  policy,
      KOKKOS_LAMBDA(const team_policy::member_type &team) 
    {

        const int c = team.league_rank();   // chain id


        for (long long step = 0; step < n_iters ; ++step)
        {
            // Only one thread per team chooses the move type
            Kokkos::single(Kokkos::PerTeam(team), [&](){
                auto r = pool.get_state();
                flip_data_local.flipMoveType(c) = r.drand(0., 1.);
                pool.free_state(r);
            });
            team.team_barrier();

            if (flip_data_local.flipMoveType(c) < 0.5f) {
                hierarchicalFlipMoveAddEnd(team, flip_data_local, c, pool);
            } else {
                hierarchicalFlipMoveAddStart(team, flip_data_local, c, pool);
            }
            team.team_barrier();

            const int accepted = flip_data_local.accept_move(c);  // uniform read by all threads
            // Evaluate and accept
            if (accepted) {
                hierarchicalDeltaE_1(team, flip_data_local, c);  // writes d.d_E_1(c)
            }
                team.team_barrier();
            if (accepted) {
                if (flip_data_local.flipMoveType(c) < 0.5f) {
                    hierarchicalOneKernel_AddEnd_FirstPart(team, flip_data_local, c, pool);
                } else {
                    hierarchicalOneKernel_AddStart_FirstPart(team, flip_data_local, c, pool);
                }
               
            }
            team.team_barrier();
 
        }

        // Final bookkeeping per chain
        hierarchicalEnergy(team, flip_data_local, c);               // writes d.newE(c)
        team.team_barrier();
        hierarchicalOneKernel_Reconnect(team, flip_data_local, c, pool);
        team.team_barrier();
        flip_data_local.E(c) = flip_data_local.newE(c);

    }); // end parallel_for
}

 // simpler use python 
 // too lazy rewrite here  
void XY_SAW_LongInteraction::gyration() {
}
void XY_SAW_LongInteraction::updateData() {

    int c = 0 ; 

    auto start_host = Kokkos::create_mirror_view(flip_data.start_conformation);
    Kokkos::deep_copy(start_host, flip_data.start_conformation);
    start_conformation = start_host(0);

    auto end_host = Kokkos::create_mirror_view(flip_data.end_conformation);
    Kokkos::deep_copy(end_host, flip_data.end_conformation);
    end_conformation = end_host(0);


    float r2 = lattice->radius(start_conformation, end_conformation);
    e2e_distance_2 << r2;

    auto E_host = Kokkos::create_mirror_view(flip_data.E);
    Kokkos::deep_copy(E_host, flip_data.E);
    E = E_host(0);

    energy << E;
    energy_2 << E * E;
    energy_4 << E * E * E * E;

    float sum_sin_1 = 0.0;
    float sum_cos_1 = 0.0;
    long long current = start_conformation;

    Kokkos::deep_copy(h_sequence_on_lattice_h, flip_data.sequence_on_lattice);
    Kokkos::deep_copy(h_lattice_nodes_positions_h, flip_data.lattice_nodes_positions);

    for (int e = 0; e < L; e++) {
        sum_sin_1 += sin(h_sequence_on_lattice_h(c, h_lattice_nodes_positions_h(c, e)));
        sum_cos_1 += cos(h_sequence_on_lattice_h(c, h_lattice_nodes_positions_h(c, e)));

    }

    sum_sin_1 /= L;
    sum_cos_1 /= L;

    mags_sin << sum_sin_1;
    mags_cos << sum_cos_1;


    magnetization_1 <<  std::sqrt(sum_sin_1 * sum_sin_1 + sum_cos_1 * sum_cos_1);
    magnetization_2 << sum_sin_1 * sum_sin_1 + sum_cos_1 * sum_cos_1;
    magnetization_4
            << (sum_sin_1 * sum_sin_1 + sum_cos_1 * sum_cos_1) * (sum_sin_1 * sum_sin_1 + sum_cos_1 * sum_cos_1);


    gyration();
    //defect();
}


void XY_SAW_LongInteraction::out_angle_data(std::fstream &out, long long n_steps) {
// Called after update; sequence is already got on host 


    for (int c =0 ; c < N_CHAINS; c++) {
        out << n_steps << " ";

        for (int e = 0; e < L; e++) {
            out <<  h_sequence_on_lattice_h(c, h_lattice_nodes_positions_h(c, e)) << " " ;
        }
        out  << std::endl;   
    }

}

void XY_SAW_LongInteraction::out_dir_data(std::fstream &out, long long n_steps) {
    // Called after update; sequence is already got on host and positions 
 
    auto ind_start = Kokkos::create_mirror_view(flip_data.start_index_in_nodes_position);
    Kokkos::deep_copy(ind_start, flip_data.start_index_in_nodes_position);
    auto E_host = Kokkos::create_mirror_view(flip_data.E);
    Kokkos::deep_copy(E_host, flip_data.E);
    //E = E_host(0);


    for (int c = 0; c < N_CHAINS; c++ ) { 
        out << n_steps << " " <<  ind_start(c) << " " << E_host(c) << " ";
        auto ls =  lattice->lattice_size();
        for (int i = 0; i < L; i++) {
            int pos = h_lattice_nodes_positions_h(c, i);
            int x = pos % ls;
            int y = (pos % (ls * ls )) / ls;
            int z = pos / ( ls * ls );
            out << x << " " << y << " " << z << " "; 
        }
        out  << std::endl;    
    }

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
    out << magnetization_1.mean() << " " << magnetization_1.errorbar() << " ";


    out << eigen1.mean() << " " << eigen1.errorbar() << " ";
    out << eigen2.mean() << " " << eigen2.errorbar() << " ";
    out << eigen3.mean() << " " << eigen3.errorbar() << " ";


    out << gyration_2_trace.mean() << " " << gyration_2_trace.errorbar() << " ";
    out << gyration_2_direct.mean() << " " << gyration_2_direct.errorbar() << " ";

    out << asphericity_collect.mean() << " " << asphericity_collect.errorbar() << " ";

    out << std::endl;

    
}