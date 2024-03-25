#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
#include <cstdlib>
#include <random>
#include <list>
#include <unordered_map>
#include <iostream>
using namespace std;

void findNeighbours(vector<vector<int>> &cell_matrix,
                    vector<vector<int>> &neighbours_matrix,
                    vector<double> &pos_array,
                    vector<double> &L_d,
                    size_t nb_moving_part, 
                    int Nx, int Ny, int Nz,
                    double h, int kappa);

void naiveAlgo(vector<vector<int>> &neighbours_matrix,
               vector<double> &pos_array,
               size_t nb_moving_part,
               double h, 
               int kappa);

void printNeighbours(vector<vector<int>> &neighbours_matrix_linked, 
                     vector<vector<int>> &neighbours_matrix_naive);

void CompareNeighbours(const std::vector<std::vector<int>> &neighbours_matrix_linked,
                     const std::vector<std::vector<int>> &neighbours_matrix_naive);
