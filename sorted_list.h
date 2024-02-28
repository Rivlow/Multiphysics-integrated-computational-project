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


void linkedListAlgo(vector<double> &particle_x, vector<double> &particle_y, vector<double> &particle_z, vector<unsigned> &particle_i, vector<unsigned> &particle_j, vector<unsigned> &particle_k, 
                 vector<vector<int>> &neighbours_matrix, double L[3], const int &Nx, const int &Ny, const int &Nz, const double &h, const int &kappa);

    

void naiveAlgo(vector<double> &particle_x, vector<double> &particle_y, vector<double> &particle_z, 
                 vector<vector<int>> &neighbours_matrix, const int &h, const int &kappa);

