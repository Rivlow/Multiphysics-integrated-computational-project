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


void findNeighbours(vector<double> &pos_arr, vector<vector<unsigned>> &cell_matrix, vector<vector<unsigned>> &neighbours_matrix, 
                    vector<double> &L_d, const unsigned &Nx, const unsigned &Ny, const unsigned &Nz, const double &h, const int &kappa);
    
void naiveAlgo(vector<double> &pos_arr, vector<vector<unsigned>> &neighbours_matrix, const double &h, const int &kappa);

void printNeighbours(vector<vector<unsigned>> &neighbours_matrix_1, vector<vector<unsigned>> &neighbours_matrix_2);

