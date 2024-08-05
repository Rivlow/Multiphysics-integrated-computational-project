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
                vector<vector<double>> &gradW,
                vector<vector<double>> &W,
                vector<vector<double>> &viscosity,
                vector<double> &nb_neighbours,
                vector<double> &type,
                vector<double> &pos);

void naiveAlgo(GeomData &geomParams,
               SimulationData &simParams, 
               vector<int> &neighbours,
               vector<double> &pos);

void printNeighbours(vector<int> &neighbours_linked,
                     vector<int> &neighbours_naive,
                     vector<double> &pos);

void CompareNeighbours( std::vector<std::vector<int>> &neighbours_matrix_linked,
                        std::vector<std::vector<int>> &neighbours_matrix_naive);

void compareAlgo(GeomData &geomParams,
                 SimulationData &simParams, 
                 vector<vector<int>> &cell_matrix,
                 vector<int> &neighbours,
                 vector<vector<double>> &gradW,
                 vector<vector<double>> &W,
                 vector<vector<double>> &viscosity,
                 vector<double> &nb_neighbours,
                 vector<double> &type,
                 vector<double> &pos);

#endif // FIND_NEIGHBOURS_H