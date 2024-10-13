//
// Created by Kamilla Faizullina on 08.04.2024.
//
#include "Model.h"
#include<iostream>
#include <random>
#include <fstream>
#include <chrono>
#include<queue>

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
    next_monomers.resize(lattice->NumberOfNodes(),NO_SAW_NODE );
    previous_monomers.resize(lattice->NumberOfNodes(), NO_SAW_NODE );
    directions.resize(lattice->NumberOfNodes(),NO_SAW_NODE ); //directions enumerated from o to dim2()
    lattice_nodes_positions.resize(number_of_spins(),NO_SAW_NODE);
}

XY_SAW_LongInteraction::XY_SAW_LongInteraction(long length) : SAW_model<double>(length) {
#ifdef REGIME_2D
    lattice = new Lattice_2D(2*L+OUT_Length);
#else
    lattice = new Lattice_3D(2*L+OUT_Length);
#endif
    if (lattice!= nullptr) {
        LatticeInitialization();
        SequenceOnLatticeInitialization();
        StartConfiguration();
    }
}

void XY_SAW_LongInteraction::SequenceOnLatticeInitialization(){
    sequence_on_lattice.resize(lattice->NumberOfNodes(),NO_XY_SPIN);
    used_coords.resize(lattice->NumberOfNodes(), false  );
}

void XY_SAW_LongInteraction::StartConfiguration() {
#ifdef STARTDEFAULT
    start_conformation=0;
    end_conformation=L-1;
    lattice_nodes_positions[0] = start_conformation;
    lattice_nodes_positions[number_of_spins()-1] = end_conformation;
    for (int i = 1; i < L-1; i++)
    {
        previous_monomers[i]=i-1;
        sequence_on_lattice[i]=PI;
        next_monomers[i]=i+1;
        lattice_nodes_positions[i] = i;
    }
    sequence_on_lattice[0] = PI;
    sequence_on_lattice[end_conformation] = PI; //начальная последовательность
    next_monomers[0] = 1;
    previous_monomers[L-1] = L-2;
    //E =  -(L-1); //Hamiltonian out of (n-1) pairs of spins
    for (int i = 0; i < L-1; i++)
    {
        directions[i]=0; //all directions are the right moves
    }
#elif STARTHALF
    coord_t middle = number_of_spins()/2 - 1;
    start_conformation=0;
    lattice_nodes_positions[0] = start_conformation;
    //First part
    int i_pos = 0;
    for (int i = 1; i < middle; i++)
    {
        previous_monomers[i]=lattice->map_of_contacts_int[lattice->ndim2()*i +1];
        sequence_on_lattice[i]=PI;
        next_monomers[i]=lattice->map_of_contacts_int[lattice->ndim2()*i +0];
        directions[i]=0;
        lattice_nodes_positions[i] = i;
    }
    //middle
    i_pos = middle;
    previous_monomers[middle]=lattice->map_of_contacts_int[lattice->ndim2()*middle + 1];
    sequence_on_lattice[middle]=PI;
    next_monomers[middle]=lattice->map_of_contacts_int[lattice->ndim2()*middle + 2];
    directions[middle] = 2; //Go Up
    lattice_nodes_positions[i_pos] = middle;
    i_pos += 1;
    middle = next_monomers[middle];
    previous_monomers[middle]=lattice->map_of_contacts_int[lattice->ndim2()*(middle)+ 3];
    sequence_on_lattice[middle]=PI;
    next_monomers[middle]=lattice->map_of_contacts_int[lattice->ndim2()*(middle) + 1];
    directions[middle] = 1;
    lattice_nodes_positions[i_pos] = middle;
    i_pos += 1;
    middle = next_monomers[middle];
    lattice_nodes_positions[i_pos] = middle;
    for (int i = number_of_spins()/2 + 2; i < number_of_spins() ; i++)
    {
        previous_monomers[middle ]=lattice->map_of_contacts_int[lattice->ndim2()*middle  +0];
        sequence_on_lattice[middle ]=PI;
        next_monomers[middle ]=lattice->map_of_contacts_int[lattice->ndim2()*middle  +1];
        directions[middle] = 1;
        middle = next_monomers[middle];
        lattice_nodes_positions[i] = middle;
    }
    end_conformation=middle;
    sequence_on_lattice[0] = PI;
    sequence_on_lattice[end_conformation] = PI; //начальная последовательность
    next_monomers[0] = lattice->map_of_contacts_int[lattice->ndim2()*0 +0];;
    previous_monomers[end_conformation] = lattice->map_of_contacts_int[lattice->ndim2()*end_conformation +0]; ;
    lattice_nodes_positions[number_of_spins() - 1] = end_conformation;
    std::fstream myStream;
    std::string filename = "For_Debug_AStart.out";
    myStream.open(filename,std::fstream::out);
    for (long j =0; j < number_of_spins() ; j++){
       myStream << lattice_nodes_positions[j] << " ";
    }
    myStream << std::endl;
    for (long j =0; j< lattice->NumberOfNodes() ; j++){
        myStream << j << " " << next_monomers[j] << " " << previous_monomers[j] << " " <<  sequence_on_lattice[j];
        myStream << std::endl;
    }
    myStream << std::endl;
    myStream.close();
#endif
    E = Energy();
}

double XY_SAW_LongInteraction::Energy() {
    double H = 0;
    double r;

    for ( long i = 0; i < number_of_spins(); i++) {
        for ( long j = i + 1; j < number_of_spins(); j++) {
            r = lattice->radius(lattice_nodes_positions[i], lattice_nodes_positions[j]);
            r = std::pow(r, R_POWER/2.0);
            H = H + (std::cos(sequence_on_lattice[lattice_nodes_positions[i]]-sequence_on_lattice[lattice_nodes_positions[j]]))/ r;
        }
    }
    return -H;
}

std::uniform_real_distribution<double> distribution_urd(0.0,1.0);
#ifdef SEED
std::mt19937 generator(URD_SEED+1);
#else
std::mt19937 generator(std::chrono::steady_clock::now().time_since_epoch().count());
#endif

void XY_SAW_LongInteraction::FlipMove_AddEnd(long direction, double spinValue) {

    coord_t new_point = lattice->map_of_contacts_int[lattice->ndim2() * end_conformation + direction];
    double oldspin = sequence_on_lattice[start_conformation];

    //self-avoidance condition:
    if (sequence_on_lattice[new_point] != NO_XY_SPIN) return;

    coord_t save_start_conformation;

    // delete the beginning of SAW
    save_start_conformation = start_conformation;
    start_conformation = next_monomers[start_conformation];
    next_monomers[save_start_conformation] = NO_SAW_NODE;
    previous_monomers[start_conformation] = NO_SAW_NODE;
    sequence_on_lattice[save_start_conformation] = NO_XY_SPIN;

    //add the new monomer at the end of SAW
    next_monomers[end_conformation] = new_point;
    sequence_on_lattice[new_point] = spinValue; //new spin value
    previous_monomers[new_point] = end_conformation;
    end_conformation = new_point;

    for (int i = 1; i < number_of_spins(); i++) {
        lattice_nodes_positions[i-1] = lattice_nodes_positions[i];
    }
    lattice_nodes_positions[number_of_spins()-1] = end_conformation;
    double new_E = Energy();

    double p1 = exp( -( J * (new_E - E) ));
    double p_metropolis = std::min(1.0, p1);

    double q_ifaccept = distribution_urd(generator) ;

    if (q_ifaccept < p_metropolis) { // accept the new state
        E = new_E;
        sequence_on_lattice[save_start_conformation] = NO_XY_SPIN;
        directions[save_start_conformation] = NO_SAW_NODE;
        directions[previous_monomers[end_conformation]] = direction;
    }
    else {
        //reject new state
        //delete end
        coord_t del = end_conformation;
        end_conformation = previous_monomers[end_conformation];
        next_monomers[end_conformation] = NO_SAW_NODE;
        previous_monomers[del] = NO_SAW_NODE;
        sequence_on_lattice[del] = NO_XY_SPIN;

        //add the previous beginning
        previous_monomers[start_conformation] = save_start_conformation;
        next_monomers[save_start_conformation] = start_conformation;
        start_conformation = save_start_conformation;
        sequence_on_lattice[start_conformation] = oldspin;

        for (int i = number_of_spins() - 1 ; i > 0; i--) {
            lattice_nodes_positions[i] = lattice_nodes_positions[i-1];
        }
        lattice_nodes_positions[0] = start_conformation;
    }
}

void XY_SAW_LongInteraction::FlipMove_AddStart(long direction, double spinValue) {

    coord_t new_point = lattice->map_of_contacts_int[lattice->ndim2() * start_conformation + direction];
    double oldspin = sequence_on_lattice[end_conformation];

    //self-avoidance condition:
    if (sequence_on_lattice[new_point] != NO_XY_SPIN) return;

    coord_t save_end_conformation;

    //delete end
    save_end_conformation = end_conformation;
    end_conformation = previous_monomers[end_conformation];
    previous_monomers[save_end_conformation] = NO_SAW_NODE;
    next_monomers[end_conformation] = NO_SAW_NODE;
    sequence_on_lattice[save_end_conformation] = NO_XY_SPIN;

    //add the new beginning
    previous_monomers[start_conformation] = new_point;
    sequence_on_lattice[new_point] = spinValue; //выбор спина
    next_monomers[new_point] = start_conformation;
    start_conformation = new_point;

    for (int i = number_of_spins() - 1 ; i > 0; i--) {
        lattice_nodes_positions[i] = lattice_nodes_positions[i-1];
    }
    lattice_nodes_positions[0] = start_conformation;
    double new_E = Energy();

    double p1 = exp( -( J * (new_E - E) ));
    double p_metropolis = std::min(1.0, p1);
    double q_ifaccept = distribution_urd(generator) ;

    if (q_ifaccept < p_metropolis) {
        E = new_E;
        sequence_on_lattice[save_end_conformation] = NO_XY_SPIN;
        directions[end_conformation] = NO_SAW_NODE;
        directions[start_conformation] = lattice->inverse_steps[direction];
    }
    else {
        //reject the new state
        //delete starte
        coord_t del = start_conformation;
        start_conformation = next_monomers[start_conformation];
        previous_monomers[start_conformation] = NO_SAW_NODE;
        next_monomers[del] = NO_SAW_NODE;
        sequence_on_lattice[del] = NO_XY_SPIN;

        //readd the end of the saw
        next_monomers[end_conformation] = save_end_conformation;
        previous_monomers[save_end_conformation] = end_conformation;
        end_conformation = save_end_conformation;
        sequence_on_lattice[end_conformation] = oldspin;

        for (int i = 1; i < number_of_spins(); i++) {
            lattice_nodes_positions[i-1] = lattice_nodes_positions[i];
        }
        lattice_nodes_positions[number_of_spins()-1] = end_conformation;
    }
}

void XY_SAW_LongInteraction::ClusterStep(double flipdirection) {
    std::uniform_int_distribution<long int> distribution_spin(0, L-1);
#ifdef SEED
    std::mt19937 generator_spin(URD_SEED+10);
#else
    std::mt19937 generator_spin(std::chrono::steady_clock::now().time_since_epoch().count());
#endif
    //long choose_spin = distribution_spin(generator_spin);
}

template<>
void SAW_model<double>::Reconnect(short direction) {
    long c = 0; //
    coord_t step_coord = lattice->map_of_contacts_int[lattice->ndim2()*end_conformation + direction];

    // test self avoidance condition
    if (sequence_on_lattice[step_coord] == NO_XY_SPIN ||
        next_monomers[step_coord] == NO_SAW_NODE ||
        step_coord == previous_monomers[end_conformation])
    {
        return;
    }

    long new_end = next_monomers[step_coord];
    next_monomers[step_coord]=end_conformation;
    directions[step_coord]=lattice->inverse_steps[direction];
    c = end_conformation;
    long int new_c;
    while (c!=new_end)
    {
        new_c=previous_monomers[c];
        next_monomers[c]=previous_monomers[c];
        directions[c]=lattice->inverse_steps[directions[new_c]];
        c=new_c;
    }
    long int temp_prev_next = next_monomers[new_end];
    previous_monomers[end_conformation]=step_coord;
    c=end_conformation;
    while (c!=new_end)
    {
        new_c=next_monomers[c];
        previous_monomers[new_c]=c;
        c=new_c;
    }
    end_conformation=new_end;
    previous_monomers[new_end]=temp_prev_next;
    next_monomers[new_end]=NO_SAW_NODE;
    directions[new_end]=NO_SAW_NODE;

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
    energy_2 << E*E;
    energy_4 << E*E*E*E;

    double sum_sin_1 = 0.0;
    double sum_cos_1 = 0.0;
    long int current = start_conformation;
    for (int e = 0; e < L ; e++)
    {
        sum_sin_1  += sin(sequence_on_lattice[current]);
        sum_cos_1  += cos(sequence_on_lattice[current]);
        current = next_monomers[current];
    }

    sum_sin_1/=L;
    sum_cos_1/=L;

    mags_sin << sum_sin_1;
    mags_cos << sum_cos_1;

    magnetization_2 <<sum_sin_1*sum_sin_1 + sum_cos_1*sum_cos_1;
    magnetization_4 << (sum_sin_1*sum_sin_1 + sum_cos_1*sum_cos_1)*(sum_sin_1*sum_sin_1 + sum_cos_1*sum_cos_1);

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