//
// Created by Kamilla Faizullina on 08.04.2024.
//

#ifndef INTERACTION_SAW_MODELS_MODEL_H
#define INTERACTION_SAW_MODELS_MODEL_H

#include <eigen3/Eigen/Dense>

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

#define N_CHAINS 10 

//const int N_CHAINS = 10; 

const float PI = std::atan(1.0)*4;

struct FlipMoveData {
    // Device-accessible data from Lattice
    Kokkos::View<int *, Kokkos::CudaSpace> map_of_contacts_int;
    Kokkos::View<int *, Kokkos::CudaSpace> inverse_steps;
    int ndim2;

    // Device-accessible data from Model
    Kokkos::View<float **, Kokkos::CudaSpace> sequence_on_lattice;
    Kokkos::View<coord_t **, Kokkos::CudaSpace> next_monomers;
    Kokkos::View<coord_t **, Kokkos::CudaSpace> previous_monomers;
    Kokkos::View<short **, Kokkos::CudaSpace> directions;
    Kokkos::View<coord_t **, Kokkos::CudaSpace> lattice_nodes_positions;

    Kokkos::View<coord_t *, Kokkos::CudaSpace> x_coords;
    Kokkos::View<coord_t *, Kokkos::CudaSpace> y_coords;
    Kokkos::View<coord_t *, Kokkos::CudaSpace> z_coords;


    // Scalars
    //float J;
    Kokkos::View<float, Kokkos::CudaSpace> J;
    //float E;

    Kokkos::View<float*, Kokkos::CudaSpace> E;
    Kokkos::View<float*, Kokkos::CudaSpace> newE;

    Kokkos::View<coord_t*, Kokkos::CudaSpace> start_conformation;
    Kokkos::View<coord_t*, Kokkos::CudaSpace> end_conformation;

    //int end_conformation;
    //int start_conformation;
    //int L;
    //int lattice_side_device;

    Kokkos::View<int, Kokkos::CudaSpace> L;
    Kokkos::View<int, Kokkos::CudaSpace> lattice_side_device;
    // Random number generator pool
    Kokkos::Random_XorShift64_Pool <Kokkos::Cuda> rand_pool;

    Kokkos::View<float*, Kokkos::CudaSpace> oldspin;
    Kokkos::View<int*, Kokkos::CudaSpace> oldIndex; // not index --- it is really coord 
    Kokkos::View<int*, Kokkos::CudaSpace> newIndex; // not index --- it is really coord 
    Kokkos::View<coord_t*, Kokkos::CudaSpace> save_start_conformation;
    Kokkos::View<coord_t*, Kokkos::CudaSpace> save_end_conformation;

    Kokkos::View<coord_t*, Kokkos::CudaSpace> start_index_in_nodes_position;

    Kokkos::View<int*, Kokkos::CudaSpace> direction;
    Kokkos::View<float*, Kokkos::CudaSpace> spinValue;

    Kokkos::View<float, Kokkos::CudaSpace> PI;
    //PI = std::atan(1.0)*4;

    Kokkos::View<int*, Kokkos::CudaSpace> i_index;
    Kokkos::View<int*, Kokkos::CudaSpace> j_index;

    Kokkos::View<int, Kokkos::CudaSpace> N_pairs;

    Kokkos::View<bool*, Kokkos::CudaSpace> accept_move; 

    Kokkos::View<float*, Kokkos::CudaSpace> flipMoveType;

    Kokkos::View<float**> localField;  // shape: (N,2)
    Kokkos::View<float*>  fieldNorm;   // shape: (N)

    Kokkos::View<int, Kokkos::CudaSpace> spin_relax  ;
    Kokkos::View<int*, Kokkos::CudaSpace> chosenIndices_relax  ;

    Kokkos::View<float*, Kokkos::CudaSpace> d_E_1; 
    Kokkos::View<float*, Kokkos::CudaSpace> d_E_2; 

    //int n_chains = 10;
};


class Model {
public:
    //Model() {};
    //Model (int length);
    KOKKOS_INLINE_FUNCTION int number_of_spins() {return L;}
    KOKKOS_INLINE_FUNCTION short ndim2() {
        //short dim2 = -1;
        //(lattice!= nullptr) ? dim2 = lattice->ndim2() : dim2 = -1;
        return lattice->ndim2();
    }
    Lattice *lattice = nullptr;
    void set_J (float J_) {J = J_;}
//protected:
    //Model-specific Energy function; returns float as J is expected to be float also
    //KOKKOS_FUNCTION virtual void Energy () = 0;
    //KOKKOS_FUNCTION virtual float Energy_Add_Start () = 0;
    //KOKKOS_FUNCTION virtual float Energy_Add_End () = 0;

    int L; //Length of the model chain
    float E; //current value for energy; float as J
    float J; //Interaction Energy
    Kokkos::View<int*, Kokkos::CudaSpace> lattice_side; //("lattice_side", 1);
    Kokkos::View<int*, Kokkos::HostSpace>::HostMirror lattice_side_host;

    FlipMoveData flip_data;

    Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> rand_pool;




};

// Abstract Class for geometry related work
template<class SpinType>
class SAW_model : public Model {
public:
    SAW_model() {};
    SAW_model(int length);

 //   KOKKOS_INLINE_FUNCTION virtual void Reconnect(short direction) = 0; //Only Geometry changes --- the same for all SAW Models



    //KOKKOS_INLINE_FUNCTION
    //virtual void FlipMove_AddEnd () = 0; //depends on spin variables
    //virtual void FlipMove_AddStart () = 0;
    virtual void runMCMCOnDevice(long long MC_STEPS) = 0;
    //KOKKOS_INLINE_FUNCTION virtual void FlipMove_AddEnd1 () = 0; //depends on spin variables
    //KOKKOS_INLINE_FUNCTION virtual void FlipMove_AddStart1 () = 0; //depends on spin variables

    //KOKKOS_INLINE_FUNCTION virtual void FlipMove_AddEnd (int direction, SpinType spinvalue) = 0; //depends on spin variables
    //KOKKOS_INLINE_FUNCTION virtual void FlipMove_AddStart (int direction, SpinType spinvalue) = 0; //depends on spin variables
   // virtual void ClusterStep (float flipdirection) = 0; //depends on spin variables

    void LatticeInitialization();

    //virtual void (std::fstream& out) = 0 ;
    virtual void out_MC_data(std::fstream& out, long long n_steps) = 0 ;
    virtual void updateData() = 0;

//protected:
    std::valarray<SpinType> sequence_on_lattice_h;
    typename Kokkos::View<SpinType**, Kokkos::CudaSpace>::HostMirror h_sequence_on_lattice_h;
    std::valarray<int> next_monomers_h;
    Kokkos::View<int**, Kokkos::CudaSpace> next_monomers;
    Kokkos::View<int**, Kokkos::CudaSpace>  previous_monomers;
    std::valarray<int> previous_monomers_h;
    int end_conformation = 0;
    int start_conformation = 0;
    std::valarray<short> directions_h; // n-1 edges of SAW on the lattice; //directions enumerated from o to dim2()
    Kokkos::View<short**, Kokkos::CudaSpace>  directions;

    mc_stats::ScalarObservable<float> e2e_distance_2;
    mc_stats::ScalarObservable<float> gyration_2_trace;
    mc_stats::ScalarObservable<float> gyration_2_direct;

    int* lattice_nodes_positions_h;
    Kokkos::View<int**, Kokkos::CudaSpace>::HostMirror h_lattice_nodes_positions_h;

    Kokkos::View<int**, Kokkos::CudaSpace>::HostMirror h_next_monomers_h;
    Kokkos::View<int**, Kokkos::CudaSpace>::HostMirror h_previous_monomers_h;
    Kokkos::View<short**, Kokkos::CudaSpace>::HostMirror h_directions_h;

    Kokkos::View<int**, Kokkos::CudaSpace> lattice_nodes_positions;
    Kokkos::View<SpinType**, Kokkos::CudaSpace> sequence_on_lattice;

    //static Kokkos::Random_XorShift64_Pool<Kokkos::Cuda> rand_pool_host;

    //void initializePool(unsigned int seed=17);
};

//Class for XY int-interacting Model on SAWs
class XY_SAW_LongInteraction : public  SAW_model<float> {
public:
    XY_SAW_LongInteraction() {};
    XY_SAW_LongInteraction(int length, float J);

    //KOKKOS_INLINE_FUNCTION void Reconnect(short direction); //Only Geometry changes --- the same for all SAW Models



    //KOKKOS_INLINE_FUNCTION void FlipMove_AddEnd1 () override;
    //KOKKOS_INLINE_FUNCTION
    //void FlipMove_AddEnd () override;
    //void FlipMove_AddStart () override;
    void runMCMCOnDevice(long long n_steps) override;
    //KOKKOS_INLINE_FUNCTION void FlipMove_AddStart1() override;

//    KOKKOS_INLINE_FUNCTION void FlipMove_AddEnd (int direction, float spinValue) override;
//    KOKKOS_INLINE_FUNCTION void FlipMove_AddStart(int direction, float spinValue) override;
   // KOKKOS_INLINE_FUNCTION void ClusterStep (float flipdirection);

    void SequenceOnLatticeInitialization();
    void StartConfiguration();


    void gyration(); 
   // void defect(std::fstream &out, long long n_steps);
    void out_MC_data(std::fstream& out, long long n_steps);
    void updateData();
    void out_angle_data(std::fstream& out, long long n_steps);
    void out_dir_data(std::fstream& out, long long n_steps);

//protected:
    std::valarray<bool> used_coords;

   // KOKKOS_FUNCTION void Energy ();
    //KOKKOS_FUNCTION float Energy_Add_Start () ;
    //KOKKOS_FUNCTION float Energy_Add_End () ;


    mc_stats::ScalarObservable<float> energy;
    mc_stats::ScalarObservable<float> energy_2;
    mc_stats::ScalarObservable<float> energy_4;

    mc_stats::ScalarObservable<float> mags_sin;
    mc_stats::ScalarObservable<float> mags_cos;
    mc_stats::ScalarObservable<float> magnetization_1;
    mc_stats::ScalarObservable<float> magnetization_2;
    mc_stats::ScalarObservable<float> magnetization_4;

    mc_stats::ScalarObservable<float> eigen1;
    mc_stats::ScalarObservable<float> eigen2;
    mc_stats::ScalarObservable<float> eigen3;

    mc_stats::ScalarObservable<float> asphericity_collect;
 



};


struct Vector3 {
    float x, y, z;
  
    KOKKOS_INLINE_FUNCTION
    Vector3() : x(0.0), y(0.0), z(0.0) {}
  
    // Overload the += operator for accumulation.
    KOKKOS_INLINE_FUNCTION
    Vector3& operator+=(const Vector3& rhs) {
      x += rhs.x;
      y += rhs.y;
      z += rhs.z;
      return *this;
    }
  };

// Structure to hold the gyration tensor




// Custom struct to accumulate the independent components of a symmetric 3x3 gyration tensor.
struct GyrationTensor {
    float q00, q01, q02, q11, q12, q22;
  
    KOKKOS_INLINE_FUNCTION
    GyrationTensor() : q00(0.0), q01(0.0), q02(0.0),
                       q11(0.0), q12(0.0), q22(0.0) {}
  
    // Overload the += operator for reduction.
    KOKKOS_INLINE_FUNCTION
    GyrationTensor& operator+=(const GyrationTensor& rhs) {
      q00 += rhs.q00;
      q01 += rhs.q01;
      q02 += rhs.q02;
      q11 += rhs.q11;
      q12 += rhs.q12;
      q22 += rhs.q22;
      return *this;
    }
  };
  
  // Specialize Kokkos::reduction_identity for GyrationTensor so that it can be used with Kokkos::Sum.
  /*namespace Kokkos {
  template <>
  struct reduction_identity<GyrationTensor> {
    KOKKOS_INLINE_FUNCTION
    static GyrationTensor max() { return GyrationTensor(); }
  };
  }*/ // namespace Kokkos
  


  namespace Kokkos {
    template <>
    struct reduction_identity<Vector3> {
      KOKKOS_INLINE_FUNCTION
      static Vector3 sum() { return Vector3(); }
    };
  
    template <>
    struct reduction_identity<GyrationTensor> {
      KOKKOS_INLINE_FUNCTION
      static GyrationTensor sum() { return GyrationTensor(); }
    };
  }


  KOKKOS_INLINE_FUNCTION
float angleDiff(float a, float b) {
  float diff = a - b;
  while(diff > M_PI)  diff -= 2.0 * M_PI;
  while(diff < -M_PI) diff += 2.0 * M_PI;
  return diff;
}


#endif //INTERACTION_SAW_MODELS_MODEL_H
