//
// Created by Kamilla Faizullina on 15.04.2024.
//

#include "ExactEnumeration.h"
#include <random>
#include "Model.h"

using namespace std;

vector<vector<tuple<int, int>>> get_all_conformations(int length){
    static vector<tuple<int,int>> steps = {make_tuple(1, 0), make_tuple(-1, 0), make_tuple(0, 1),  make_tuple(0, -1)};
    vector<vector<tuple<int, int>>> result;
    vector<tuple<int, int>> temp;
    tuple<int, int> new_point;
    if(length==1){
        result.push_back({make_tuple(0, 0)});
        return result;
    }
    else{
        result=get_all_conformations(length-1);
        vector<vector<tuple<int, int>>>  new_conformations;
        for (int i=0; i<result.size(); i++){
            for (tuple<int, int> step : steps){
                //new_point = make_tuple(get<0>(result[i][-1])+get<0>(step), get<1>(result[i][-1])+get<1>(step) );
                new_point = make_tuple(get<0>(result[i].back())+get<0>(step), get<1>(result[i].back())+get<1>(step) );
                if(find(result[i].begin(), result[i].end(), new_point)!=result[i].end() ){
                    continue;
                }
                temp = result[i];
                temp.push_back(new_point);
                new_conformations.push_back(temp);
                temp.clear();
            }
        }
        return new_conformations;
    }
}


double radius2(const tuple<int, int>& start, const tuple<int, int>& end) {
    int xdiff = get<0>(end) - get<0>(start);
    int ydiff = get<1>(end) - get<1>(start);
    return xdiff*xdiff + ydiff*ydiff;
}

double Hamiltonian(const vector<tuple<int, int>>&  conformation,
                   const vector<double>& spin_sequence,
                   double J) {
    double r;
    double H = 0.;
    int i = 0;
    //int check_number = 0;
    for (auto out : conformation) {
        int j = 0;
        for (auto inner : conformation) {
            //if (out!=inner) {
            if (i!=j) {
                r = radius2(out, inner);
                r = std::sqrt(r);
                //if (r < 1.001) {
#ifdef XYSHORT
                if (r < 1.001)
#endif
                    H = H + (cos(spin_sequence[i] - spin_sequence[j])) / (r * r * r);
                //}
            }
            j+=1;
        }
        i+=1;
    }
    return -1.*H/2.;
}


void MeanValues(int L) {
    long double Z = 0;
    long double E = 0;
    long double R = 0;
    long double r2;

    vector<vector<tuple<int, int>>> conformations = get_all_conformations(L);
    vector<double> spin_sequence;
    std::uniform_real_distribution<double> distribution_theta(0, 2.0 * PI);
    std::mt19937 generators_theta;
#ifdef SEED
    generators_theta.seed(103);
#else
    generators_theta.seed(chrono::steady_clock::now().time_since_epoch().count());
#endif
    for (int spin_n = 0; spin_n < L; spin_n++) {
        spin_sequence.push_back(distribution_theta(generators_theta));
    }
    double weight;
    std::cout << "L J E r2" << std::endl;

    for (int Attempt = 0 ; Attempt < 10; Attempt++) {
        for (double J = 0.0; J <= 1.5; J += 0.3) {
            Z = 0;
            R = 0;
            E = 0;
            for (vector<tuple<int, int>> conformation: conformations) {
                r2 = radius2(conformation.front(),
                             conformation.back());
                for (int n_conf = 0; n_conf < 10000; n_conf++) {

                    double H = Hamiltonian(conformation, spin_sequence, J);
                    weight = exp(-J * H);
                    Z = Z + weight;
                    R = R + r2 * weight;
                    E = E + H * weight;

                    // Get new spin sequence
                    for (int spin_n = 0; spin_n < L; spin_n++) {
                        spin_sequence[spin_n] = distribution_theta(generators_theta);
                    }
                }
            }
            std::cout << L << " " << J << " " << E / Z << " " << R / Z << std::endl;
        }
    }

}