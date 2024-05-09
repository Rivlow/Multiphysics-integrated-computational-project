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

void sortedList(GeomData &geomParams,
                SimulationData &simParams, 
                vector<vector<int>> &cell_matrix,
                vector<int> &neighbours,
                vector<vector<double>> &gradW_matrix,
                vector<vector<double>> &artificial_visc_matrix,
                vector<double> &nb_neighbours,
                vector<double> &type,
                vector<double> &pos);

void naiveAlgo(GeomData &geomParams,
               SimulationData &simParams, 
               vector<vector<int>> &neighbours_matrix,
               vector<double> &pos);

void printNeighbours(vector<vector<int>> &neighbours_matrix_linked, 
                     vector<vector<int>> &neighbours_matrix_naive);

void CompareNeighbours( std::vector<std::vector<int>> &neighbours_matrix_linked,
                        std::vector<std::vector<int>> &neighbours_matrix_naive);

#endif // FIND_NEIGHBOURS_H