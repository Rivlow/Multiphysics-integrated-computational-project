#ifndef CUBE_H
#define CUBE_H
#include <vector>

using namespace std;

int evaluateNumberParticles(vector<double> &L, double &s);

void meshcube(vector<double> &o, vector<double> &L, double &s, std::vector<double> &pos_arr);

void clearAllVectors(vector<vector<double>> & artificial_visc_matrix, vector<vector<unsigned>> &neighbours_matrix, 
                     vector<vector<unsigned>> &cell_matrix, vector<vector<double>> &gradW_matrix);

void meshBoundary(vector<double> &o, vector<double> &L, double &s, std::vector<double> &bound_arr);

#endif // CUBE_H