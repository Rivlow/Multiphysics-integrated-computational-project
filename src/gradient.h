#ifndef GRADIENT_H
#define GRADIENT_H

#include <stdio.h>
#include <vector>
#include "structure.h"


using namespace std;

void gradW(const SimulationData& params, 
           vector<vector<double>> &gradW_matrix,
           vector<vector<int>> &neighbours_matrix,
           vector<double> &pos_array,
           int nb_moving_part,
           int Nx, int Ny, int Nz);

void setSpeedOfSound(const SimulationData& params,
                     vector<double> &c_array,
                     vector<double> &rho_array);

void setPressure(const SimulationData& params,
                 vector<double> &p_array,
                 vector<double> &rho_array, 
                 int nb_moving_part);

void setArtificialViscosity(const SimulationData& params,
                            int t,
                            vector<vector<double>> &artificial_visc_matrix,
                            vector<vector<int>> &neighbours_matrix,
                            vector<double> &c_array,
                            vector<double> &pos_array,
                            vector<double> &rho_array,
                            vector<double> &u_array, 
                            int nb_moving_part);

void continuityEquation(const SimulationData& params,
                        vector<vector<int>> &neighbours_matrix,
                        vector<vector<double>> &gradW_matrix,
                        vector<double> &pos_array,
                        vector<double> &u_array,
                        vector<double> &drhodt_array,
                        vector<double> &rho_array,
                        vector<double> &mass_array, 
                        int nb_moving_part);

void momentumEquation(const SimulationData& params,
                      int t,
                      vector<vector<int>> &neighbours_matrix,
                      vector<vector<double>> &gradW_matrix,
                      vector<vector<double>> &artificial_visc_matrix,
                      vector<double> &mass_array,
                      vector<double> &dudt_array,
                      vector<double> &rho_array,
                      vector<double> &p_array, 
                      vector<double> &c_array,
                      vector<double> &pos_array,
                      vector<double> &u_array,
                      int nb_moving_part);





void update(const SimulationData& params,
            vector<double> &pos_array,
            vector<double> &u_array,
            vector<double> &rho_array,
            vector<double> &drhodt_array,
            vector<double> &dudt_array);

#endif // GRADIENT_H