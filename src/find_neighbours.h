#ifndef FIND_NEIGHBOURS_H
#define FIND_NEIGHBOURS_H

#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
#include <cstdlib>
#include <random>
#include <list>
#include <unordered_map>
#include <iostream>

#include "structure.h"

using namespace std;

void sorted_list(const SimulationData& params, 
                 vector<vector<int>> &cell_matrix,
                 vector<vector<int>> &neighbours_matrix,
                 vector<double> &pos_array,
                 int nb_moving_part, 
                 int Nx, int Ny, int Nz);

void naiveAlgo(const SimulationData& params, 
               vector<vector<int>> &neighbours_matrix,
               vector<double> &pos_array,
               int nb_moving_part);

void printNeighbours(vector<vector<int>> &neighbours_matrix_linked, 
                     vector<vector<int>> &neighbours_matrix_naive);

void CompareNeighbours(const std::vector<std::vector<int>> &neighbours_matrix_linked,
                       const std::vector<std::vector<int>> &neighbours_matrix_naive);

#endif // FIND_NEIGHBOURS_H