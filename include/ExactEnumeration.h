//
// Created by Kamilla Faizullina on 15.04.2024.
//

#ifndef INTERACTION_SAW_MODELS_EXACTENUMERATION_H
#define INTERACTION_SAW_MODELS_EXACTENUMERATION_H

#include<vector>
#include<map>
#include <tuple>
#include <iostream>
#include <algorithm>
#include <fstream>

#ifdef REGIME_2D
std::vector <int> vector_for_distance(std::vector<std::tuple<int, int>> saw );
std::vector<std::vector<std::tuple<int, int>>> get_all_conformations(int length);
#else
std::vector <int> vector_for_distance(std::vector<std::tuple<int, int, int>> saw );
std::vector<std::vector<std::tuple<int, int, int>>> get_all_conformations(int length);
#endif

void MeanValues(int L);

#endif //INTERACTION_SAW_MODELS_EXACTENUMERATION_H
