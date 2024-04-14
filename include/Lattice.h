#ifndef INTERACTION_SAW_MODELS_LATTICE_H
#define INTERACTION_SAW_MODELS_LATTICE_H

#include <valarray>
#include <vector>

typedef unsigned long coord_t;

class Lattice {
public:
    //Lattice();
    Lattice(long max_seq_size = 10) { lattice_side = max_seq_size; };

    inline const virtual short int ndim() = 0 ;
    inline const virtual short int ndim2() = 0;

    inline long int lattice_size() {return lattice_side;};
    inline long int NumberOfNodes () {return number_of_nodes;};

    std::valarray<coord_t> map_of_contacts_int;
    std::valarray<int> inverse_steps;
protected:
    virtual void create_lattice() = 0;

    long lattice_side = 0;
    int number_of_nodes = 0;
    std::vector<std::vector<int> > steps;
};


class Lattice_2D : public Lattice {
public:
    Lattice_2D(long max_seq_size = 10);
    inline const short int ndim() { return 2; }
    inline const short int ndim2() {return 4;}

private:
    void create_lattice();
};


class Lattice_3D : public Lattice {
public:
    Lattice_3D(long max_seq_size = 10);
    inline const short int ndim() { return 3; }
    inline const short int ndim2() {return 6;}

private:
    void create_lattice();
};



#endif //INTERACTION_SAW_MODELS_LATTICE_H
