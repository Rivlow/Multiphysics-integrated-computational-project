#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <stdio.h>
#include <vector>
#include <string>

#include "tools.h"
#include "structure.h"


using namespace std;

void initializeMass( SimulationData& params, 
                    vector<double> &rho_array,
                    vector<double> &mass_array);

void initializeRho( SimulationData& params,
                   vector<double> &pos_array,
                   vector<double> &rho_array);

void initializeVelocity( SimulationData& params, 
                        vector<double> &u_array);

void initializeViscosity( SimulationData& params, 
                         vector<vector<double>> &artificial_visc_matrix);

void checkTimeStep(SimulationData &params, 
                   int t,
                   vector<double> pos,
                   vector<double> c,
                   vector<vector<int>> &neighbours_matrix,
                   vector<vector<double>> &artificial_visc_matrix);

#endif // INITIALIZE_H