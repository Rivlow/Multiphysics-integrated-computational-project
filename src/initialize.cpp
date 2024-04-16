#include <stdio.h>
#include <vector>
#include <string>
#include <limits>
#include <iostream>

#include "initialize.h"
#include "structure.h"
#include "gradient.h"



using namespace std;

void initializeMass(GeomData &geomParams,
                    SimulationData &simParams, 
                    vector<double> &rho,
                    vector<double> &mass){

    double s = geomParams.s;  
    bool PRINT = simParams.PRINT;  
    double V = s * s * s;
    int nb_fixed_part = simParams.nb_fixed_part;

    for (int i = 0; i < nb_fixed_part; i++){
        mass[i] = rho[i] * V;
    }

    if (PRINT){
        cout << "initializeMass passed" << endl;
    }
}

void initializeRho(ThermoData &thermoParams,
                   SimulationData &simParams,
                   vector<double> &pos,
                   vector<double> &rho){

    string state_initial_condition = simParams.state_initial_condition;
    string state_equation = simParams.state_equation;
    bool PRINT = simParams.PRINT;
    int nb_moving_part = simParams.nb_moving_part;

    double R = thermoParams.R;
    double T = thermoParams.T;
    double M = thermoParams.M;
    double rho_0 = thermoParams.rho_0;
    double rho_fixed = thermoParams.rho_fixed;
    double rho_moving = thermoParams.rho_moving;
    double c_0 = thermoParams.c_0;
    double gamma = thermoParams.gamma;
    double g = thermoParams.g;

    int nb_fixed_part = simParams.nb_fixed_part;


    if (state_initial_condition == "Hydrostatic"){
        if (state_equation == "Ideal gaz law"){
            for (int i = 0; i < nb_fixed_part; i++){

                rho[i] = (i < nb_moving_part) ? rho_0 * (1 + M * rho_0 
                                * g * pos[3 * i + 2] / (R * T)) : rho_fixed;
            }
        }

        if (state_equation == "Quasi incompresible fluid"){

            double B = c_0 * c_0 * rho_0 / gamma;
            for (int i = 0; i < nb_fixed_part; i++){

                rho[i] = (i < nb_moving_part) ? rho_0 * (1 + rho_0 
                                * g * pos[3 * i + 2] / B) : rho_fixed;
            }
        }
    }
    else{

        for (int i = 0; i < nb_fixed_part; i++){
            rho[i] = (i < nb_moving_part) ? rho_moving : rho_fixed;
        }
    }

    if (PRINT){

        cout << "initializeRho passed" << endl;
    }
}

void initializeVelocity(ThermoData &thermoParams,
                        SimulationData &simParams, 
                        vector<double> &u){


    bool PRINT = simParams.PRINT;
    int nb_moving_part = simParams.nb_moving_part;

    for (int i = 0; i < nb_moving_part; i++){

        u[3 * i] = simParams.u_init[0];
        u[3 * i + 1] = simParams.u_init[1];
        u[3 * i + 2] = simParams.u_init[2];
    }

    if (PRINT){
       cout << "initializeVelocity passed" << endl;
    }
}

void initializeViscosity(SimulationData &simParams, 
                         vector<vector<double>> &pi_matrix){

    bool PRINT = simParams.PRINT;
    int size_pi_matrix = pi_matrix.size();

    for (int i = 0; i < size_pi_matrix; i++){

        int size_artificial_visc = pi_matrix[i].size();

        for (int j = 0; j < size_artificial_visc; j++){
            pi_matrix[i][j] = 0.0;
        }
    }

    if (PRINT){
        cout << "initializeViscosity passed" << endl;
    }
}

void checkTimeStep(GeomData &geomParams,    
                   ThermoData &thermoParams,
                   SimulationData &simParams, 
                   int t,
                   vector<double> pos,
                   vector<double> c,
                   vector<vector<int>> &neighbours_matrix,
                   vector<double> &nb_neighbours,
                   vector<vector<double>> &pi_matrix){

    double alpha = thermoParams.alpha;
    double beta = thermoParams.beta;
    double h = geomParams.h;
    double g = thermoParams.g;
    int nb_moving_part = simParams.nb_moving_part;

    double dt_f = h / abs(g);
    double dt_cv;
    double min_a = numeric_limits<double>::max();
    double max_b = numeric_limits<double>::min();

    if (t == 0){

        double prev_dt = simParams.dt;
        double dt_f = h / abs(g);
        simParams.dt = (simParams.dt > dt_f) ? dt_f : simParams.dt;
        double next_dt = simParams.dt;

        if (abs(prev_dt - next_dt) != 0){
            cout << "dt has to be modified, was : " << prev_dt << " and is now : " << next_dt << endl;
        }
    }
    else{

        for (int n = 0; n < nb_moving_part; n++){

            vector<double> &artificial_visc = pi_matrix[n];
            vector<int> &neighbours = neighbours_matrix[n];

            double c_a = c[n];
            int size_neighbours = neighbours.size();

            for (int idx = 0; idx < size_neighbours; idx++){

                double pi_ab = artificial_visc[idx];
                max_b = (pi_ab > max_b) ? pi_ab: max_b;

            }

            double val = h/(c_a + 0.6*(alpha*c_a + beta*max_b));
            min_a = (val < min_a) ? val : min_a;

        }

        dt_cv = min_a;
        double dt_final = min(0.25*dt_f, 0.4*dt_cv);

        string state_equation = simParams.state_equation;

        if (state_equation == "Ideal gaz law"){

            double prev_dt = simParams.dt;
            simParams.dt = (dt_final < simParams.dt) ? dt_final : simParams.dt;
            double next_dt = simParams.dt;
            
            if (abs(prev_dt - next_dt) != 0){
                cout << "dt has to be modified, was " << prev_dt << " and is now " << next_dt << endl;
            }
        }
        else{

            double prev_dt = simParams.dt;
            simParams.dt = (dt_final < simParams.dt) ? dt_final : simParams.dt;
            double next_dt = simParams.dt;

            if (abs(prev_dt - next_dt) != 0){
                cout << "dt has to be modified, was " << prev_dt << " and is now " << next_dt << endl;
            }
        }
    }
}
