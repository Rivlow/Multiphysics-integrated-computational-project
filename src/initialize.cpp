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

    double s = params.s;  
    bool PRINT = params.PRINT;  
    double V = s * s * s;

    for (int i = 0; i < int(rho_array.size()); i++)
    {
        mass_array[i] = rho_array[i] * V;
    }

    if (PRINT){
        cout << "initializeMass passed" << endl;
    }
}

void initializeRho(const SimulationData& params,
                   vector<double> &pos_array,
                   vector<double> &rho_array){

    string state_initial_condition = params.state_initial_condition;
    string state_equation = params.state_equation;
    double R = params.R;
    double T = params.T;
    double M = params.M;
    double rho_0 = params.rho_0;
    double rho_fixed = params.rho_fixed;
    double rho_moving = params.rho_moving;
    double c_0 = params.c_0;
    double gamma = params.gamma;
    double g = params.g;
    bool PRINT = params.PRINT;
    int nb_moving_part = params.nb_moving_part;



    if (state_initial_condition == "Hydrostatic"){
        if (state_equation == "Ideal gaz law"){
            for (int i = 0; i < int(rho_array.size()); i++){
                rho_array[i] = (i < nb_moving_part) ? rho_0 * (1 + M * rho_0 
                                * g * pos_array[3 * i + 2] / (R * T)) : rho_fixed;
            }
        }
        if (params.state_equation == "Quasi incompresible fluid"){

            double B = c_0 * c_0 * rho_0 / gamma;
            for (int i = 0; i < int(rho_array.size()); i++){
                rho_array[i] = (i < nb_moving_part) ? rho_0 * (1 + rho_0 
                                * g * pos_array[3 * i + 2] / B) : rho_fixed;
            }
        }
    }
    else{

        for (int i = 0; i < int(rho_array.size()); i++){
            rho_array[i] = (i < nb_moving_part) ? rho_moving : rho_fixed;
        }
    }

    if (PRINT){

        cout << "initializeRho passed" << endl;
    }
}

void initializeVelocity(const SimulationData& params, 
                        vector<double> &u_array){


    bool PRINT = params.PRINT;
    int nb_moving_part = params.nb_moving_part;

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

    if (PRINT){
       cout << "initializeVelocity passed" << endl;
    }
}

void initializeViscosity(const SimulationData& params, 
                         vector<vector<double>> &artificial_visc_matrix){

    bool PRINT = params.PRINT;

    for (int i = 0; i < int(artificial_visc_matrix.size()); i++){
        for (int j = 0; j < int(artificial_visc_matrix[i].size()); j++){
            artificial_visc_matrix[i][j] = 0.0;
        }
    }

    if (PRINT){
        cout << "initializeViscosity passed" << endl;
    }
}
