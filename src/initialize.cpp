#include <stdio.h>
#include <vector>
#include <string>
#include "initialize.h"
#include <iostream>

using namespace std;

void initializeMass(const SimulationData& params, vector<double> &rho_array,
                    vector<double> &mass_array){
    double V = params.s * params.s * params.s;
    for (size_t i = 0; i < rho_array.size(); i++)
    {
        mass_array[i] = rho_array[i] * V;
    }

    if (params.PRINT){
        cout << "initializeMass passed" << endl;
    }
}

void initializeRho(const SimulationData& params,
                   vector<double> &pos_array,
                   vector<double> &rho_array, 
                   size_t nb_moving_part){

    if (params.state_initial_condition == "Hydrostatic")
    {
        if (params.state_equation == "Ideal gaz law")
        {
            for (size_t i = 0; i < rho_array.size(); i++)
            {
                rho_array[i] = (i < nb_moving_part) ? params.rho_0 * (1 + params.M * params.rho_0 * params.g * pos_array[3 * i + 2] / (params.R * params.T)) : params.rho_fixed;
            }
        }
        if (params.state_equation == "Quasi incompresible fluid")
        {

            double B = params.c_0 * params.c_0 * params.rho_0 / params.gamma;
            for (size_t i = 0; i < rho_array.size(); i++)
            {
                rho_array[i] = (i < nb_moving_part) ? params.rho_0 * (1 + params.rho_0 * params.g * pos_array[3 * i + 2] / B) : params.rho_fixed;
            }
        }
    }
    else
    {
        for (size_t i = 0; i < rho_array.size(); i++)
        {
            rho_array[i] = (i < nb_moving_part) ? params.rho_moving : params.rho_fixed;
        }
    }

    if (params.PRINT){
        
        cout << "initializeRho passed" << endl;
    }
}

void initializeVelocity(const SimulationData& params, vector<double> &u_array,
                        vector<double> &u_init,
                        size_t nb_moving_part){

    for (size_t i = 0; i < u_array.size() / 3; i++){

        if (i < nb_moving_part){
            u_array[3 * i] = u_init[0];
            u_array[3 * i + 1] = u_init[1];
            u_array[3 * i + 2] = u_init[2];
        }
        else{
            u_array[3 * i] = 0;
            u_array[3 * i + 1] = 0;
            u_array[3 * i + 2] = 0;
        }
    }

    if (params.PRINT){
       cout << "initializeVelocity passed" << endl;
    }
}

void initializeViscosity(const SimulationData& params, vector<vector<double>> &artificial_visc_matrix)
{
    for (size_t i = 0; i < artificial_visc_matrix.size(); i++){
        for (size_t j = 0; j < artificial_visc_matrix[i].size(); j++)
        {
            artificial_visc_matrix[i][j] = 0.0;
        }
    }

    if (params.PRINT){
        cout << "initializeViscosity passed" << endl;
    }
}
