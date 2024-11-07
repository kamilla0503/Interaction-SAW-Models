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

XY_SAW_LongInteraction::XY_SAW_LongInteraction(long length) : SAW_model<double>(length) {
#ifdef REGIME_2D
    lattice = new Lattice_2D(2 * L + OUT_Length);
#else
    lattice = new Lattice_3D(2*L+OUT_Length);
#endif
    if (lattice != nullptr) {
        LatticeInitialization();
        SequenceOnLatticeInitialization();
        StartConfiguration();
    }
    rand_pool = Kokkos::Random_XorShift64_Pool<Kokkos::DefaultExecutionSpace>(/*seed=*/12345);
    Kokkos::fence();
    printf("Finish all configs \n");

    //lattice_side = lattice->lattice_side;
    //rand_pool.init(12345,256);
}

void XY_SAW_LongInteraction::SequenceOnLatticeInitialization() {
    //sequence_on_lattice_h = new double [lattice->NumberOfNodes()]{NO_XY_SPIN};
    sequence_on_lattice_h.resize(lattice->NumberOfNodes(), NO_XY_SPIN);
    used_coords.resize(lattice->NumberOfNodes(), false);
}

void XY_SAW_LongInteraction::StartConfiguration() {
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
        sequence_on_lattice_h[i] = PI;
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
        sequence_on_lattice_h[i]=PI;
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
        sequence_on_lattice_h[middle ]=PI;
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

    E = Energy();

    std::cout << "Model creation after energy" << std::endl;
    printf("Energy after all = %f \n", E);

// Copy Kokkos::View members from Lattice
    flip_data.map_of_contacts_int = lattice->map_of_contacts_int;
    flip_data.inverse_steps = lattice->inverse_steps;
    flip_data.ndim2 = lattice->ndim2();
    printf("Finish lattice fields = %f \n", E);
// Copy Kokkos::View members from XY_SAW_LongInteraction
    flip_data.sequence_on_lattice = sequence_on_lattice;
    printf("Finish  seq = %f \n", E);
    flip_data.next_monomers = next_monomers;
    printf("Finish next = %f \n", E);
    flip_data.previous_monomers = previous_monomers;
    printf("Finish prevs = %f \n", E);
    flip_data.directions = directions;
    printf("Finish directions = %f \n", E);
    flip_data.lattice_nodes_positions = lattice_nodes_positions;
    printf("Finish nodes positions = %f \n", E);

// Scalars
    //flip_data.NO_SAW_NODE = NO_SAW_NODE;
    //flip_data.NO_XY_SPIN = NO_XY_SPIN;
    flip_data.J = J;
    printf("Finish J = %f \n", E);
    flip_data.E = E;
    printf("Finish E= %f \n", E);

// Random pool
    Kokkos::Random_XorShift64_Pool <Kokkos::Cuda> rand_pool = Kokkos::Random_XorShift64_Pool<Kokkos::Cuda>(/* seed or execution space */);
    flip_data.rand_pool_ptr = &rand_pool;

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

// Define the Energy function using Kokkos parallelization
KOKKOS_FUNCTION
double XY_SAW_LongInteraction::Energy() {
    double H = 0.0;  // Total energy
    auto local_L = L;
    auto lattice_side_local = lattice_side_host(0);
    auto lattice_nodes_positions_local = lattice_nodes_positions;
    auto sequence_on_lattice_local = sequence_on_lattice;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Kokkos::Cuda>(0, local_L), KOKKOS_LAMBDA(
    const long i,
    double &local_H) {
        double r;
        double energy_i = 0.0;  // Local energy contribution for this i
        if (i < local_L ) {
            for (long j = i + 1; j < local_L; j++) {
                r = radius(lattice_nodes_positions_local(i), lattice_nodes_positions_local(j),
                           lattice_side_local);
                r = Kokkos::pow(r, R_POWER / 2.0);
                energy_i += Kokkos::cos(
                        sequence_on_lattice_local(lattice_nodes_positions_local(i)) -
                        sequence_on_lattice_local(lattice_nodes_positions_local(j))) /
                                r;
            }
        }
        local_H += energy_i;  // Add local energy contribution to the reduction variable
    }, H);  // H is the total energy accumulated across all threads
    return -H;  // Return negative of the total energy
}

std::uniform_real_distribution<double> distribution_urd(0.0, 1.0);
#ifdef SEED
std::mt19937 generator(URD_SEED + 1);
#else
std::mt19937 generator(std::chrono::steady_clock::now().time_since_epoch().count());
#endif

KOKKOS_INLINE_FUNCTION
void XY_SAW_LongInteraction::FlipMove_AddEnd(long direction, double spinValue) {
    printf("FlipMove_AddEnd \n");
    coord_t new_point = lattice->map_of_contacts_int(lattice->ndim2() * end_conformation + direction);
    //coord_t new_point = flip_data.map_of_contacts_int(lattice->ndim2() * end_conformation + direction);
    //std::cout << "new_point " << new_point << std::endl;
    printf("FlipMove_AddEnd new point = %ld \n", new_point);
    double oldspin = sequence_on_lattice(start_conformation);
    //std::cout << "oldspin " << oldspin << std::endl;
    printf("FlipMove_AddEnd new point = %f \n", oldspin);
    //self-avoidance condition:
    if (sequence_on_lattice(new_point) != NO_XY_SPIN) return;
    //std::cout << "sequence_on_lattice(new_point)" << sequence_on_lattice(new_point) << std::endl;
    coord_t save_start_conformation;

    // delete the beginning of SAW
    save_start_conformation = start_conformation;
    printf("FlipMove_AddEnd save start = %ld \n", save_start_conformation);
    start_conformation = next_monomers(start_conformation);
    printf("FlipMove_AddEnd start = %ld \n", start_conformation);
    next_monomers(save_start_conformation) = NO_SAW_NODE;
    previous_monomers(start_conformation) = NO_SAW_NODE;
    sequence_on_lattice(save_start_conformation) = NO_XY_SPIN;

    //add the new monomer at the end of SAW
    next_monomers(end_conformation) = new_point;
    sequence_on_lattice(new_point) = spinValue; //new spin value
    previous_monomers(new_point) = end_conformation;
    end_conformation = new_point;
    //std::cout << "Movement of positions  " << start_conformation << " " << spinValue << std::endl;

    for (int i = 1; i < number_of_spins(); i++) {
        lattice_nodes_positions(i - 1) = lattice_nodes_positions(i);
    }
    lattice_nodes_positions(this->number_of_spins() - 1) = end_conformation;

    double new_E = Energy();

    double p1 = exp(-(J * (new_E - E)));
    double p_metropolis = Kokkos::min(1.0, p1);

    auto rand_gen = rand_pool.get_state();
    // Generate a random number between 0.0 and 1.0
    double q_ifaccept = rand_gen.drand(0., 1.);
   if (q_ifaccept < p_metropolis) { // accept the new state
        E = new_E;
        sequence_on_lattice(save_start_conformation) = NO_XY_SPIN;
        directions(save_start_conformation) = NO_SAW_NODE;
        directions(previous_monomers(end_conformation)) = direction;
    } else {
        //reject new state
        //delete end
        coord_t del = end_conformation;
        end_conformation = previous_monomers(end_conformation);
        next_monomers(end_conformation) = NO_SAW_NODE;
        previous_monomers(del) = NO_SAW_NODE;
        sequence_on_lattice(del) = NO_XY_SPIN;

        //add the previous beginning
        previous_monomers(start_conformation) = save_start_conformation;
        next_monomers(save_start_conformation) = start_conformation;
        start_conformation = save_start_conformation;
        sequence_on_lattice(start_conformation) = oldspin;

        for (int i = this->number_of_spins() - 1; i > 0; i--) {
            lattice_nodes_positions(i) = lattice_nodes_positions(i - 1);
        }
        lattice_nodes_positions(0) = start_conformation;
    }
    rand_pool.free_state(rand_gen);
}

KOKKOS_INLINE_FUNCTION
void XY_SAW_LongInteraction::FlipMove_AddStart(long direction, double spinValue) {

    coord_t new_point = lattice->map_of_contacts_int(lattice->ndim2() * start_conformation + direction);
    double oldspin = sequence_on_lattice(end_conformation);

    if (sequence_on_lattice(new_point) != NO_XY_SPIN) return;

    coord_t save_end_conformation;

    //delete end
    save_end_conformation = end_conformation;
    end_conformation = previous_monomers(end_conformation);
    previous_monomers(save_end_conformation) = NO_SAW_NODE;
    next_monomers(end_conformation) = NO_SAW_NODE;
    sequence_on_lattice(save_end_conformation) = NO_XY_SPIN;

    //add the new beginning
    previous_monomers(start_conformation) = new_point;
    sequence_on_lattice(new_point) = spinValue; //выбор спина
    next_monomers(new_point) = start_conformation;
    start_conformation = new_point;

    for (int i = this->number_of_spins() - 1; i > 0; i--) {
        lattice_nodes_positions(i) = lattice_nodes_positions(i - 1);
    }
    lattice_nodes_positions(0) = start_conformation;
    double new_E = Energy();

    double p1 = exp(-(J * (new_E - E)));
    double p_metropolis = Kokkos::min(1.0, p1);

    auto rand_gen = rand_pool.get_state();
    double q_ifaccept = rand_gen.drand(0., 1.);

    if (q_ifaccept < p_metropolis) {
        E = new_E;
        sequence_on_lattice(save_end_conformation) = NO_XY_SPIN;
        directions(end_conformation) = NO_SAW_NODE;
        directions(start_conformation) = lattice->inverse_steps(direction);
    } else {
        //reject the new state
        //delete starte
        coord_t del = start_conformation;
        start_conformation = next_monomers(start_conformation);
        previous_monomers(start_conformation) = NO_SAW_NODE;
        next_monomers(del) = NO_SAW_NODE;
        sequence_on_lattice(del) = NO_XY_SPIN;

        //readd the end of the saw
        next_monomers(end_conformation) = save_end_conformation;
        previous_monomers(save_end_conformation) = end_conformation;
        end_conformation = save_end_conformation;
        sequence_on_lattice(end_conformation) = oldspin;

        for (int i = 1; i < this->number_of_spins(); i++) {
            lattice_nodes_positions(i - 1) = lattice_nodes_positions(i);
        }
        lattice_nodes_positions(this->number_of_spins() - 1) = end_conformation;
    }

    rand_pool.free_state(rand_gen);

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
    double r2 = lattice->radius(start_conformation, end_conformation);
    e2e_distance_2 << r2;
    energy << E;
    energy_2 << E * E;
    energy_4 << E * E * E * E;

    double sum_sin_1 = 0.0;
    double sum_cos_1 = 0.0;
    long int current = start_conformation;
    for (int e = 0; e < L; e++) {
        sum_sin_1 += sin(sequence_on_lattice[current]);
        sum_cos_1 += cos(sequence_on_lattice[current]);
        current = next_monomers[current];
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