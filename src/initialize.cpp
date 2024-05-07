#include <stdio.h>
#include <vector>
#include <string>
#include <limits>
#include <iostream>

#include "initialize.h"
#include "structure.h"
#include "gradient.h"
#include "tools.h"

using namespace std;


void initMass(GeomData &geomParams,
              SimulationData &simParams, 
              vector<double> &rho,
              vector<double> &mass){

    double s = geomParams.s;  
    bool PRINT = simParams.PRINT;  
    double V = s * s * s;
    int rho_size = rho.size();

    #pragma omp parallel for   
    for (int i = 0; i < rho_size; i++)
        mass[i] = rho[i] * V;
    

    if (PRINT) cout << "initMass passed" << endl;
    
}

void initRho(ThermoData &thermoParams,
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

    int rho_size = rho.size();


    if (state_initial_condition == "Hydrostatic"){
        if (state_equation == "Ideal gaz law"){

            #pragma omp parallel for   
            for (int i = 0; i < rho_size; i++)
                rho[i] = (i < nb_moving_part) ? rho_0 * (1 + M * rho_0 
                                * g * pos[3 * i + 2] / (R * T)) : rho_fixed;    
        }

        if (state_equation == "Quasi incompresible fluid"){

            double B = c_0 * c_0 * rho_0 / gamma;
            for (int i = 0; i < rho_size; i++)
                rho[i] = (i < nb_moving_part) ? rho_0 * (1 + rho_0 
                                * g * pos[3 * i + 2] / B) : rho_fixed;  
        }
    }
    else{

        #pragma omp parallel for   
        for (int i = 0; i < rho_size; i++)
            rho[i] = (i < nb_moving_part) ? rho_moving : rho_fixed;
        
    }

    if (PRINT) cout << "initRho passed" << endl;
    
}

void initVelocity(ThermoData &thermoParams,
                        SimulationData &simParams, 
                        vector<double> &u){


    bool PRINT = simParams.PRINT;
    int nb_moving_part = simParams.nb_moving_part;
    int u_size = u.size();

    #pragma omp parallel for   
    for (int i = 0; i < u_size / 3; i++){

        if (i < nb_moving_part){
            u[3 * i] = simParams.u_init[0];
            u[3 * i + 1] = simParams.u_init[1];
            u[3 * i + 2] = simParams.u_init[2];
        }
        else{
            u[3 * i] = 0;
            u[3 * i + 1] = 0;
            u[3 * i + 2] = 0;
        }
    }

    if (PRINT){
       cout << "initVelocity passed" << endl;
    }
}

void initViscosity(SimulationData &simParams, 
                         vector<vector<double>> &pi_matrix){

    bool PRINT = simParams.PRINT;
    int size_pi_matrix = pi_matrix.size();
    #pragma omp parallel for   
    for (int i = 0; i < size_pi_matrix; i++){

        int size_artificial_visc = pi_matrix[i].size();

        for (int j = 0; j < size_artificial_visc; j++){
            pi_matrix[i][j] = 0.0;
        }
    }

    if (PRINT){
        cout << "initViscosity passed" << endl;
    }
}

void checkTimeStep(GeomData &geomParams,    
                   ThermoData &thermoParams,
                   SimulationData &simParams, 
                   vector<double> &pos,
                   vector<double> &u,
                   vector<double> &c,
                   vector<vector<int>> &neighbours_matrix,
                   vector<double> &nb_neighbours,
                   vector<vector<double>> &pi_matrix){

    double alpha = thermoParams.alpha;
    double beta = thermoParams.beta;
    double h = geomParams.h;
    double g = (simParams.is_gravity) ? -9.81 : 0;
    int nb_moving_part = simParams.nb_moving_part;
    int t = simParams.t;

    double dt_f = h / abs(g);
    double dt_cv;
    double min_a = numeric_limits<double>::max();
    double max_b = numeric_limits<double>::min();

    if (t == 0){

        double prev_dt = simParams.dt;
        double F_st_max = thermoParams.F_st_max;
        double dt_f = h / sqrt(F_st_max*F_st_max + g*g);
        simParams.dt = (simParams.dt > dt_f) ? dt_f : simParams.dt;
        double next_dt = simParams.dt;

        if (abs(prev_dt - next_dt) != 0){
            cout << "dt has to be modified (timestep :" << t <<")"<<", was : " << prev_dt << " and is now : " << next_dt << endl;
        }
    }
    else{

        #pragma omp parallel for   
        for (int n = 0; n < nb_moving_part; n++){

            vector<int> &neighbours = neighbours_matrix[n];
            int size_neighbours = nb_neighbours[n];

            // Iteration over each associated neighbours
            for (int idx = 0; idx < size_neighbours; idx++){

                int i_neig = neighbours[idx];
                vector<double> rel_displ(3), rel_vel(3);

                rel_displ[0] = (pos[3 * n + 0] - pos[3 * i_neig + 0]);
                rel_displ[1] = (pos[3 * n + 1] - pos[3 * i_neig + 1]);
                rel_displ[2] = (pos[3 * n + 2] - pos[3 * i_neig + 2]);

                rel_vel[0] = (u[3 * n + 0] - u[3 * i_neig + 0]);
                rel_vel[1] = (u[3 * n + 1] - u[3 * i_neig + 1]);
                rel_vel[2] = (u[3 * n + 2] - u[3 * i_neig + 2]);

                double u_ab_x_ab = 0, x_ab_2 = 0;

                // Dot product
                for (int cord = 0; cord < 3; cord++){
                    u_ab_x_ab += rel_vel[cord] * rel_displ[cord];
                    x_ab_2 += rel_displ[cord] * rel_displ[cord];
                }

                double nu_2 = 0.01 * h * h;
                double mu_ab = (h * u_ab_x_ab) / (x_ab_2 + nu_2);
                max_b = (mu_ab > max_b)? mu_ab : max_b;
            }

            double c_a = c[n];
            double val = h/(c_a + 0.6*(alpha*c_a + beta*max_b));
            if (c_a < 0) cout << "c_a : " << c_a << " at timestep : " << t << endl;
            if (c_a < 0) printArray(c, c.size(), "c (in init)");

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
                cout << "dt modified (t :" << t <<")"<<", was : " << prev_dt << " and is now : " << next_dt << endl;
            }
        }
        else{

            double prev_dt = simParams.dt;
            simParams.dt = (dt_final < simParams.dt) ? dt_final : simParams.dt;
            double next_dt = simParams.dt;

            if (abs(prev_dt - next_dt) != 0){
                cout << "dt has to be modified (timestep :" << t <<")"<<", was : " << prev_dt << " and is now : " << next_dt << endl;
            }
        }
    }
}