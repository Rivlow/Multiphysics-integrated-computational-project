#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>

#include "initialize.h"
#include "structure.h"



using namespace std;

void initializeMass(const SimulationData& params, 
                    vector<double> &rho,
                    vector<double> &mass){

    double s = params.s;  
    bool PRINT = params.PRINT;  
    double V = s * s * s;
    int rho_size = rho.size();

    for (int i = 0; i < rho_size; i++){
        mass[i] = rho[i] * V;
    }

    if (PRINT){
        cout << "initializeMass passed" << endl;
    }
}

void initializeRho(const SimulationData& params,
                   vector<double> &pos,
                   vector<double> &rho){

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
    int rho_size = rho.size();


    if (state_initial_condition == "Hydrostatic"){
        if (state_equation == "Ideal gaz law"){
            for (int i = 0; i < rho_size; i++){

                rho[i] = (i < nb_moving_part) ? rho_0 * (1 + M * rho_0 
                                * g * pos[3 * i + 2] / (R * T)) : rho_fixed;
            }
        }

        if (params.state_equation == "Quasi incompresible fluid"){

            double B = c_0 * c_0 * rho_0 / gamma;
            for (int i = 0; i < rho_size; i++){

                rho[i] = (i < nb_moving_part) ? rho_0 * (1 + rho_0 
                                * g * pos[3 * i + 2] / B) : rho_fixed;
            }
        }
    }
    else{

        for (int i = 0; i < rho_size; i++){
            rho[i] = (i < nb_moving_part) ? rho_moving : rho_fixed;
        }
    }

    if (PRINT){

        cout << "initializeRho passed" << endl;
    }
}

void initializeVelocity(const SimulationData& params, 
                        vector<double> &u){


    bool PRINT = params.PRINT;
    int nb_moving_part = params.nb_moving_part;
    int u_size = u.size();

    for (int i = 0; i < u_size / 3; i++){

        if (i < nb_moving_part){
            u[3 * i] = params.u_init[0];
            u[3 * i + 1] = params.u_init[1];
            u[3 * i + 2] = params.u_init[2];
        }
        else{
            u[3 * i] = 0;
            u[3 * i + 1] = 0;
            u[3 * i + 2] = 0;
        }
    }

    if (PRINT){
       cout << "initializeVelocity passed" << endl;
    }
}

void initializeViscosity(const SimulationData& params, 
                         vector<vector<double>> &artificial_visc_matrix){

    bool PRINT = params.PRINT;
    int size_artificial_visc_matrix = artificial_visc_matrix.size();

    for (int i = 0; i < size_artificial_visc_matrix; i++){

        int size_artificial_visc = artificial_visc_matrix[i].size();

        for (int j = 0; j < size_artificial_visc; j++){
            artificial_visc_matrix[i][j] = 0.0;
        }
    }

    if (PRINT){
        cout << "initializeViscosity passed" << endl;
    }
}
