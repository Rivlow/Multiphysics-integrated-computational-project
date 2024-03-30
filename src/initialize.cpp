#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>

#include "initialize.h"
#include "structure.h"



using namespace std;

void initializeMass(const SimulationData& params, 
                    vector<double> &rho_array,
                    vector<double> &mass_array){
    double V = params.s * params.s * params.s;
    for (int i = 0; i < int(rho_array.size()); i++)
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
                   int nb_moving_part){

    if (params.state_initial_condition == "Hydrostatic")
    {
        if (params.state_equation == "Ideal gaz law")
        {
            for (int i = 0; i < int(rho_array.size()); i++)
            {
                rho_array[i] = (i < nb_moving_part) ? params.rho_0 * (1 + params.M * params.rho_0 * params.g * pos_array[3 * i + 2] / (params.R * params.T)) : params.rho_fixed;
            }
        }
        if (params.state_equation == "Quasi incompresible fluid")
        {

            double B = params.c_0 * params.c_0 * params.rho_0 / params.gamma;
            for (int i = 0; i < int(rho_array.size()); i++)
            {
                rho_array[i] = (i < nb_moving_part) ? params.rho_0 * (1 + params.rho_0 * params.g * pos_array[3 * i + 2] / B) : params.rho_fixed;
            }
        }
    }
    else
    {
        for (int i = 0; i < int(rho_array.size()); i++)
        {
            rho_array[i] = (i < nb_moving_part) ? params.rho_moving : params.rho_fixed;
        }
    }

    if (params.PRINT){
        
        cout << "initializeRho passed" << endl;
    }
}

void initializeVelocity(const SimulationData& params, 
                        vector<double> &u_array,
                        int nb_moving_part){

    for (int i = 0; i < int(u_array.size()) / 3; i++){

        if (i < nb_moving_part){
            u_array[3 * i] = params.u_init[0];
            u_array[3 * i + 1] = params.u_init[1];
            u_array[3 * i + 2] = params.u_init[2];
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

void initializeViscosity(const SimulationData& params, 
                         vector<vector<double>> &artificial_visc_matrix)
{
    for (int i = 0; i < int(artificial_visc_matrix.size()); i++){
        for (int j = 0; j < int(artificial_visc_matrix[i].size()); j++)
        {
            artificial_visc_matrix[i][j] = 0.0;
        }
    }

    if (params.PRINT){
        cout << "initializeViscosity passed" << endl;
    }
}
