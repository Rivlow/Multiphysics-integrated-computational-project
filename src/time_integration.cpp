#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
#include <cstdlib>
#include <random>
#include <list>
#include <unordered_map>
#include <iostream>
#include <chrono>
#include <algorithm>

#include "find_neighbours.h"
#include "time_integration.h"
#include "structure.h"
#include "NavierStokes.h"
#include "tools.h"
#include <omp.h>


void updateVariables(GeomData &geomParams,    
                     ThermoData &thermoParams,
                     SimulationData &simParams, 
                     vector<double> &pos,
                     vector<double> &u,
                     vector<double> &rho,
                     vector<double> &drhodt,
                     vector<double> &c,
                     vector<double> &p,
                     vector<double> &dudt,
                     vector<double> &mass,
                     vector<vector<double>> &viscosity,
                     vector<vector<double>> &gradW,
                     vector<vector<double>> &W,
                     vector<int> &neighbours,
                     vector<double> &nb_neighbours,
                     vector<double> type,
                     vector<double> &track_particle){

    bool PRINT = simParams.PRINT;
    string schemeIntegration = simParams.schemeIntegration;
    int nb_tot_part = pos.size()/3;

    if (schemeIntegration == "Euler"){

        // Compute D(rho)/Dt for all particles
        continuityEquation(simParams, neighbours, nb_neighbours, gradW, 
                           pos, u, drhodt, rho, mass);

        // Compute D(u)/Dt for moving particles
        momentumEquation(geomParams, thermoParams, simParams, neighbours, nb_neighbours, gradW, 
                         W, viscosity, mass, dudt, rho, p, c, pos, u, type, track_particle); 

        checkTimeStep(geomParams, thermoParams, simParams, pos, u, c,
                    neighbours, nb_neighbours);

        double dt = simParams.dt;

        #pragma omp parallel for   
        for (int n = 0; n < simParams.nb_tot_part; n++){

            rho[n] += dt * drhodt[n];

            if (rho[n] < 0){
                cerr << "Rho negative (" << rho[n] << ") at t:" << simParams.t << endl;
                exit(EXIT_FAILURE);
            }

            for (int coord = 0; coord < 3; coord++){

                pos[3 * n + coord] += dt * u[3 * n + coord];
                u[3 * n + coord] += dt * dudt[3 * n + coord];
            }
        }
    }

    else if (schemeIntegration == "RK22"){
    
        vector<double> u_half = u,
                       rho_half = rho,
                       pos_half = pos,
                       drhodt_half(nb_tot_part,0.0),
                       dudt_half(3*nb_tot_part,0.0);

        // First time step of RK22
       
        continuityEquation(simParams, neighbours, nb_neighbours, gradW, 
                           pos, u, drhodt, rho, mass);
        //printArray(drhodt,drhodt.size(), "drhodt");
        momentumEquation(geomParams, thermoParams, simParams, neighbours, nb_neighbours, gradW, 
                         W, viscosity, mass, dudt, rho, p, c, pos, u, type, track_particle); 
                        
        //checkTimeStep(geomParams, thermoParams, simParams, pos, u, c,
        //              neighbours, nb_neighbours);

        double theta = simParams.theta;
        double dt_half = simParams.dt/(2*theta);

        #pragma omp parallel for
        for (int n = 0; n < simParams.nb_tot_part; n++){

            rho_half[n] += dt_half * drhodt[n];
            
                
            if (rho_half[n] < 0){
                cout << "Rho negative (" << rho[n] << ") at t:" << simParams.t << endl;
                exit(1);
            }

            for (int coord = 0; coord < 3; coord++){
                
                pos_half[3 * n + coord] += dt_half * u[3 * n + coord];
                u_half[3 * n + coord] += dt_half * dudt[3 * n + coord];
            }
        }
        
        // Second time step of RK22
    
        continuityEquation(simParams, neighbours, nb_neighbours, gradW, 
                           pos_half, u_half, drhodt_half, rho_half, mass);

        momentumEquation(geomParams, thermoParams, simParams, neighbours, nb_neighbours, gradW, 
                        W, viscosity, mass, dudt_half, rho_half, p, c, pos_half, u_half, type, track_particle); 


        checkTimeStep(geomParams, thermoParams, simParams, pos, u, c,
                      neighbours, nb_neighbours);

        double dt = simParams.dt;
        #pragma omp parallel for
        for (int n = 0; n < simParams.nb_tot_part; n++){
            
            rho[n] += dt * ((1-theta)*drhodt[n] + theta*drhodt_half[n]);
            //cout << "for part " << n << " drhodt is " << drhodt[n]  << " and drhodt_half " << drhodt_half[n] <<" and rho " << rho[n]<< endl;
            for (int coord = 0; coord < 3; coord++){
                
                double u_temp = u[3 * n + coord] + dt*dudt_half[3 * n + coord];
                
                
                pos[3 * n + coord] += dt * ((1-theta)*u[3 * n + coord] + theta*u_temp);
                u[3 * n + coord] += dt *  ((1-theta)*dudt[3 * n + coord] + theta*dudt_half[3 * n + coord]);

            }
        }
        
    }
    else{
        cerr << "No scheme integration chosen" << endl;
        exit(EXIT_FAILURE);
    }

    if (PRINT) cout << "update passed" << endl;
    
}

void checkTimeStep(GeomData &geomParams,    
                   ThermoData &thermoParams,
                   SimulationData &simParams, 
                   vector<double> &pos,
                   vector<double> &u,
                   vector<double> &c,
                   vector<int> &neighbours,
                   vector<double> &nb_neighbours){

    double alpha = simParams.alpha;
    double beta = simParams.beta;
    double h = geomParams.h;
    
    int nb_moving_part = simParams.nb_moving_part;
    int t = simParams.t;

    double F_st_max = simParams.F_st_max;
    double dt_f = h / F_st_max;
    double dt_cv = 0;
    double min_a = numeric_limits<double>::max();
    double max_b = numeric_limits<double>::min();

    if (t == 0){

        double prev_dt = simParams.dt;
        simParams.dt = (simParams.dt > dt_f) ? dt_f : simParams.dt;
        double next_dt = simParams.dt;

        if (abs(prev_dt - next_dt) != 0)
            cout << "dt has to be modified (timestep :" << t <<")"<<", was : "
             << prev_dt << " and is now : " << next_dt << endl;
        
    }
    else{

        #pragma omp parallel for   
        for (int n = 0; n < nb_moving_part; n++){

            int size_neighbours = nb_neighbours[n];

            

            if (beta != 0.0){
                // Iteration over each associated neighbours
                for (int idx = 0; idx < size_neighbours; idx++){

                    int i_neig = neighbours[100*n + idx];
                    vector<double> rel_displ(3), rel_vel(3);

                    for (int coord = 0; coord < 3; coord++){
                        rel_displ[coord] = (pos[3 * n + coord] - pos[3 * i_neig + coord]);
                        rel_vel[coord] = (u[3 * n + coord] - u[3 * i_neig + coord]);
                    }

                    double u_ab_x_ab = 0, x_ab_2 = 0;

                    // Dot product
                    for (int coord = 0; coord < 3; coord++){
                        u_ab_x_ab += rel_vel[coord] * rel_displ[coord];
                        x_ab_2 += rel_displ[coord] * rel_displ[coord];
                    }

                    double nu_2 = 0.01 * h * h;
                    double mu_ab = (h * u_ab_x_ab) / (x_ab_2 + nu_2);
                    max_b = (mu_ab > max_b)? mu_ab : max_b;
                }
            }

            else
                max_b = 0;
        
            double c_a = c[n];
            double val = geomParams.h/(c_a + 0.6*(alpha*c_a + beta*max_b));
            min_a = (val < min_a) ? val : min_a;
        }

        dt_cv = min_a;
        double dt_final = min(0.4*dt_f, 0.25*dt_cv);
        string state_equation = simParams.state_equation;

        if (state_equation == "Ideal gaz law"){

            double prev_dt = simParams.dt;
            simParams.dt = (dt_final < simParams.dt) ? dt_final : simParams.dt;
            double next_dt = simParams.dt;
            
            if (simParams.PRINT){
                if (abs(prev_dt - next_dt) != 0){
                    cout << setprecision(15);
                    cout << "dt modified (t :" << t <<")"<<", was : " << prev_dt
                        << " and is now : " << next_dt << endl;
                }
            }
        }
        
        else{

            double prev_dt = simParams.dt;
            simParams.dt = (dt_final < simParams.dt) ? dt_final : simParams.dt;
            double next_dt = simParams.dt;

            if (simParams.PRINT){
                if (abs(prev_dt - next_dt) != 0){
                    cout << setprecision(15);
                    cout << "dt modified (t :" << t <<")"<<", was : " << prev_dt
                        << " and is now : " << next_dt << endl;
                }
            }
            
        }
    }
}
