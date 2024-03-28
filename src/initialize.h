#include <stdio.h>
#include <vector>
#include <string>
#include "tools.h"
using namespace std;

void initializeMass(const SimulationData& params, 
                    vector<double> &rho_array,
                    vector<double> &mass_array);

void initializeRho(const SimulationData& params,
                   vector<double> &pos_array,
                   vector<double> &rho_array, 
                   size_t nb_moving_part);

void initializeVelocity(const SimulationData& params, 
                        vector<double> &u_array,
                        vector<double> &u_init,
                        size_t nb_moving_part);

void initializeViscosity(const SimulationData& params, 
                         vector<vector<double>> &artificial_visc_matrix);