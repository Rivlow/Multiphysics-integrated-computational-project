#include <stdio.h>
#include <vector>
#include <string>
#include "initialize.h"
#include <iostream>

using namespace std;

void initializeMass(vector<double> &rho_array,
                    vector<double> &mass_array,
                    double s,
                    const bool PRINT){
    double V = s * s * s;
    for (size_t i = 0; i < rho_array.size(); i++)
    {
        mass_array[i] = rho_array[i] * V;
    }

    if (PRINT){
        cout << "initializeMass passed" << endl;
    }
}

void initializeRho(vector<double> &pos_array,
                   vector<double> &rho_array, 
                   size_t nb_moving_part,
                   double rho_moving,
                   double rho_fixed,
                   double rho_0, double c_0,
                   double M, double g, double R,
                   double T, double gamma,
                   std::string state_equation_chosen,
                   std::string state_initial_condition,
                   const bool PRINT){

    if (state_initial_condition == "Hydrostatic")
    {
        if (state_equation_chosen == "Ideal gaz law")
        {
            for (size_t i = 0; i < rho_array.size(); i++)
            {
                rho_array[i] = (i < nb_moving_part) ? rho_0 * (1 + M * rho_0 * g * pos_array[3 * i + 2] / (R * T)) : rho_fixed;
            }
        }
        if (state_equation_chosen == "Quasi incompresible fluid")
        {

            double B = c_0 * c_0 * rho_0 / gamma;
            for (size_t i = 0; i < rho_array.size(); i++)
            {
                rho_array[i] = (i < nb_moving_part) ? rho_0 * (1 + rho_0 * g * pos_array[3 * i + 2] / B) : rho_fixed;
            }
        }
    }
    else
    {
        for (size_t i = 0; i < rho_array.size(); i++)
        {
            rho_array[i] = (i < nb_moving_part) ? rho_moving : rho_fixed;
        }
    }

    if (PRINT){
        
        cout << "initializeRho passed" << endl;
    }
}

void initializeVelocity(vector<double> &u_array,
                        vector<double> &u_init,
                        size_t nb_moving_part, 
                        const bool PRINT){

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

    if (PRINT){
       cout << "initializeVelocity passed" << endl;
    }
}

void initializeViscosity(vector<vector<double>> &artificial_visc_matrix, const bool PRINT)
{
    for (size_t i = 0; i < artificial_visc_matrix.size(); i++){
        for (size_t j = 0; j < artificial_visc_matrix[i].size(); j++)
        {
            artificial_visc_matrix[i][j] = 0.0;
        }
    }

    if (PRINT){
        cout << "initializeViscosity passed" << endl;
    }
}
