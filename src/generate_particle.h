#ifndef CUBE_H
#define CUBE_H
#include <vector>

using namespace std;

int evaluateNumberParticles(vector<double> &L, 
                            double &s);

void meshcube(vector<double> &o, 
              vector<double> &L, 
              vector<double> &pos_arr,
              vector<double> &type_arr,
              double s);

void meshBoundary(vector<double> &o_d, 
                  vector<double> &L_d, 
                  vector<double> &bound_arr, 
                  vector<double> &type_arr,
                  double s);

void clearAllVectors(vector<vector<double>> &artificial_visc_matrix,
                     vector<vector<unsigned>> &neighbours_matrix,
                     vector<vector<unsigned>> &cell_matrix,
                     vector<vector<double>> &gradW_matrix);

#endif // CUBE_H