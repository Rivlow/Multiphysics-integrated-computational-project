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


void findNeighbours(size_t nb_moving_part, vector<double> &pos_arr, vector<vector<int>> &cell_matrix, vector<vector<int>> &neighbours_matrix, 
                    vector<double> &L_d, const int &Nx, const int &Ny, const int &Nz, const double &h, const int &kappa);
    
void naiveAlgo(size_t nb_moving_part, vector<double> &pos_arr, vector<vector<int>> &neighbours_matrix, const double &h, const int &kappa);

void printNeighbours(vector<vector<int>> &neighbours_matrix_1, vector<vector<int>> &neighbours_matrix_2);

