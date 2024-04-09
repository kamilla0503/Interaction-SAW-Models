//
// Created by Kamilla Faizullina on 08.04.2024.
//
#include "Model.h"

//used to increase length of SAWs for lattice side
#ifndef OUT_Length 2
#define OUT_Length 2
#endif
//used to define lattice nodes without spins 
#ifndef NO_SPIN
#define NO_SPIN -1
#endif


SAW_model::SAW_model(long length) {
    L = length;
}

//Initialize geometry
void SAW_model::LatticeInitialization() {
    next_monomers.resize(lattice->NumberOfNodes(),NO_SPIN ); //
    previous_monomers.resize(lattice->NumberOfNodes(), NO_SPIN );
    directions.resize(lattice->NumberOfNodes(),NO_SPIN ); //directions enumerated from o to dim2()
}

XY_SAW_LongInteraction::XY_SAW_LongInteraction(long length) : SAW_model(length) {
#ifdef REGIME_2D
    lattice = new Lattice_2D(L+OUT_Length);
#else
    lattice = new Lattice_3D(L+OUT_Length);
#endif
    if (lattice!= nullptr) {
        LatticeInitialization();
    }
}

void XY_SAW_LongInteraction::Energy() {

}

void XY_SAW_LongInteraction::FlipMove() {

}

void XY_SAW_LongInteraction::ClusterStep() {

}

void SAW_model::Reconnect(int j) {

}