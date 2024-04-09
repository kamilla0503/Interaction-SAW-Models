//
// Created by Kamilla Faizullina on 08.04.2024.
//
#include "Model.h"
#include<iostream>
//used to increase length of SAWs for lattice side
#ifndef OUT_Length 2
#define OUT_Length 2
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

SAW_model::SAW_model(long length) {
    L = length;
}

//Initialize geometry
void SAW_model::LatticeInitialization() {
    next_monomers.resize(lattice->NumberOfNodes(),NO_SAW_NODE );
    previous_monomers.resize(lattice->NumberOfNodes(), NO_SAW_NODE );
    directions.resize(lattice->NumberOfNodes(),NO_SAW_NODE ); //directions enumerated from o to dim2()
}

XY_SAW_LongInteraction::XY_SAW_LongInteraction(long length) : SAW_model(length) {
#ifdef REGIME_2D
    lattice = new Lattice_2D(L+OUT_Length);
#else
    lattice = new Lattice_3D(L+OUT_Length);
#endif
    if (lattice!= nullptr) {
        LatticeInitialization();
        SequenceOnLatticeInitialization();
        StartConfiguration();
    }
}

void XY_SAW_LongInteraction::SequenceOnLatticeInitialization(){
    sequence_on_lattice.resize(lattice->NumberOfNodes(),NO_XY_SPIN);
}

void XY_SAW_LongInteraction::StartConfiguration() {
#ifdef STARTDEFAULT
    start_conformation=0;
    end_conformation=L-1;
    for (int i = 1; i < L-1; i++)
    {
        previous_monomers[i]=i-1;
        sequence_on_lattice[i]=PI;
        next_monomers[i]=i+1;
    }
    sequence_on_lattice[0] = PI;
    sequence_on_lattice[end_conformation] = PI; //начальная последовательность
    next_monomers[0] = 1;
    previous_monomers[L-1] = L-2;
    E =  -(L-1); //Hamiltonian out of (n-1) pairs of spins
    for (int i = 0; i < L-1; i++)
    {
        directions[i]=0; //all directions are the right moves
    }
#else
#endif
}

void XY_SAW_LongInteraction::Energy() {

}

void XY_SAW_LongInteraction::FlipMove() {

}

void XY_SAW_LongInteraction::ClusterStep() {

}

void SAW_model::Reconnect(int j) {

}