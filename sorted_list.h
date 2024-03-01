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


void linkedListAlgo(vector<double> &position, vector<vector<unsigned>> &cell_pos, vector<vector<unsigned>> &neighbours_matrix, double L[3], const unsigned &Nx, const unsigned &Ny, const unsigned &Nz, const double &h, const int &kappa);
    
void naiveAlgo(vector<double> &position, vector<vector<unsigned>> &neighbours_matrix, const double &h, const int &kappa);

void printNeighbours(vector<vector<unsigned>> &neighbours_matrix_1, vector<vector<unsigned>> &neighbours_matrix_2)

