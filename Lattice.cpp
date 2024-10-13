//
// Created by Kamilla Faizullina on 01.04.2024.
//
#include <cuda_runtime_api.h>
#include <cuda.h>
#include<stdlib.h>
#include<fstream>
#include "Lattice.h"

Lattice_2D::Lattice_2D(long max_seq_size) : Lattice(max_seq_size) {
    number_of_nodes = lattice_side*lattice_side;
    create_lattice();
}

Lattice_3D::Lattice_3D(long max_seq_size) : Lattice(max_seq_size) {
    number_of_nodes = lattice_side*lattice_side*lattice_side;
    create_lattice();
}

void Lattice_2D::create_lattice() {
    long int x, y;
    ldiv_t n;
    //map_of_contacts_int.resize(lattice_side*lattice_side*ndim2());
    map_of_contacts_int = new long [lattice_side*lattice_side*ndim2()];
    for (long i =0; i<number_of_nodes ; i++){
        map_of_contacts_int[ndim2()*i] = i+1;
        map_of_contacts_int[ndim2()*i+1] = i-1;
        map_of_contacts_int[ndim2()*i+2] = i+lattice_side;
        map_of_contacts_int[ndim2()*i+3] = i-lattice_side;
        n=div(i, lattice_side);
        x=n.rem;
        y=n.quot;
        for (int j =0; j<ndim2(); j++){
            if(x==0){
                map_of_contacts_int[ndim2()*i+1] = i+lattice_side-1;
            }
            if(x==(lattice_side-1)){
                map_of_contacts_int[ndim2()*i] = i-(lattice_side-1);
            }
            if(y==0){
                map_of_contacts_int[ndim2()*i+3] = lattice_side*(lattice_side-1)+x;
            }
            if(y==(lattice_side-1)){
                map_of_contacts_int[ndim2()*i+2] = x;
            }
        }
    }
    inverse_steps.resize(ndim2());
    inverse_steps = {1,0,3,2};
#ifdef CHEKMAP2D
    std::fstream myStream;
    myStream.open("For_Debug_MapOfContacts.out",std::fstream::out);
    for (long i =0; i<number_of_nodes ; i++){
        myStream << i ;
        for (int j = 0; j < ndim2(); j++) {
            myStream << " " << map_of_contacts_int[i*ndim2()+j];
        }
        myStream << std::endl;
    }
#endif
}

void Lattice_3D::create_lattice() {
    coord_t x, y,z;
    ldiv_t n;
    //map_of_contacts_int.resize(lattice_side*lattice_side*lattice_side*ndim2());
    map_of_contacts_int = new long [lattice_side*lattice_side*ndim2()];
    coord_t l;
    for (long i =0; i<number_of_nodes ; i++){
        map_of_contacts_int[ndim2() * i] = i + 1;
        map_of_contacts_int[ndim2() * i + 1] = i - 1;
        map_of_contacts_int[ndim2() * i + 2] = i + lattice_side;
        map_of_contacts_int[ndim2() * i + 3] = i - lattice_side;
        map_of_contacts_int[ndim2() * i + 4] = i + lattice_side * lattice_side;
        map_of_contacts_int[ndim2() * i + 5] = i - lattice_side * lattice_side;

        l = lattice_side * lattice_side;
        n = div(i, l);
        z = n.quot;
        n = div( n.rem, lattice_side);
        x = n.rem;
        y = n.quot;

        for (int j = 0; j < ndim2(); j++) {
            if (x == 0) {
                map_of_contacts_int[6 * i + 1] = i + lattice_side - 1;
            }
            if (x == (lattice_side - 1)) {
                map_of_contacts_int[6 * i] = i - (lattice_side - 1);
            }
            if (y == 0) {
                map_of_contacts_int[6 * i + 3] =  lattice_side * (lattice_side - 1) + i;
            }
            if (y == (lattice_side - 1)) {
                map_of_contacts_int[6 * i + 2] = i - lattice_side * (lattice_side - 1) ;
            }
            if (z == 0) {
                map_of_contacts_int[6 * i + 5] = i + lattice_side * lattice_side * (lattice_side - 1);
            }
            if (z == lattice_side - 1) {
                map_of_contacts_int[6 * i + 4] = i - lattice_side * lattice_side * (lattice_side - 1);
            }
        }
    }
    inverse_steps.resize(ndim2());
    inverse_steps = { 1, 0, 3, 2, 5, 4 };
#ifdef CHEKMAP3D
    std::fstream myStream;
    myStream.open("For_Debug_MapOfContacts_3D.out",std::fstream::out);
    for (long i =0; i<number_of_nodes ; i++){
        myStream << i ;
        for (int j = 0; j < ndim2(); j++) {
            myStream << " " << map_of_contacts_int[i*ndim2()+j];
        }
        myStream << std::endl;
    }
#endif
    cudaMalloc(&d_map_of_contacts, lattice_side*lattice_side*ndim2()* sizeof(long));
    cudaMemcpy(d_map_of_contacts, map_of_contacts_int,
               lattice_side*lattice_side*ndim2()* sizeof(long), cudaMemcpyHostToDevice);

}

double Lattice_2D::radius(const coord_t& start, const coord_t& end) {
    long start_x = start % lattice_side;
    long start_y = start / lattice_side;
    long end_x = end % lattice_side;
    long end_y = end / lattice_side;

    //torus distance;
    double xdiff = abs(end_x - start_x);
    if (xdiff > (lattice_side/2))
        xdiff = lattice_side - xdiff;

    double ydiff = abs(end_y - start_y);
    if (ydiff > (lattice_side / 2))
        ydiff = lattice_side - ydiff;

    double r = xdiff *xdiff  + ydiff*ydiff ;

    return r;
}

double Lattice_3D::radius(const coord_t& start, const coord_t& end) {
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