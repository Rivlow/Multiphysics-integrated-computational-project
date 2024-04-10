#ifndef CUBE_H
#define CUBE_H
#include <vector>

using namespace std;

int evaluateNumberParticles(SimulationData &params);

void meshcube(SimulationData &params,
              vector<double> &pos_arr,
              vector<double> &type_arr,
              double s);

void meshBoundary(SimulationData &params,
                  vector<double> &bound_arr, 
                  vector<double> &type_arr,
                  double s);

void clearAllVectors(vector<vector<double>> &artificial_visc_matrix,
                     vector<vector<unsigned>> &neighbours_matrix,
                     vector<vector<unsigned>> &cell_matrix,
                     vector<vector<double>> &gradW_matrix);

#endif // CUBE_H