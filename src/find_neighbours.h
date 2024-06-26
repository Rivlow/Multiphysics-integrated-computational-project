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
                vector<double> &gradW,
                vector<double> &W,
                vector<double> &viscosity,
                vector<int> &nb_neighbours,
                vector<int> &type,
                vector<double> &pos,
                vector<int> &free_surface);

void naiveAlgo(GeomData &geomParams,
               SimulationData &simParams, 
               vector<int> &neighbours,
               vector<double> &pos);

void printNeighbours(vector<int> &neighbours_linked,
                     vector<int> &neighbours_naive,
                     vector<double> &pos);

void CompareNeighbours( std::vector<std::vector<int>> &neighbours_matrix_linked,
                        std::vector<std::vector<int>> &neighbours_matrix_naive);

#endif // FIND_NEIGHBOURS_H