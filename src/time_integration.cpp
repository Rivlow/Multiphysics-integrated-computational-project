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
                     vector<double> &colour,
                     vector<double> &R,
                     vector<double> &L,
                     vector<double> &N,
                     vector<double> &normal,
                     vector<double> &acc_vol,
                     vector<double> &track_particle,
                     vector<double> &Kappa,
                     vector<double> &dot_product){


    bool PRINT = simParams.PRINT;
    string scheme_integration = simParams.scheme_integration;
    int nb_tot_part = pos.size()/3;

    if (scheme_integration == "Euler"){

        // Compute D(rho)/Dt for all particles
        continuityEquation(simParams, neighbours, nb_neighbours, gradW, 
                           pos, u, drhodt, rho, mass); 

        // Compute D(u)/Dt for moving particles
        momentumEquation(geomParams, thermoParams, simParams, neighbours, nb_neighbours, gradW, 
                         W, viscosity, mass, dudt, rho, p, c, pos, u, type, colour, R, L, N, normal, acc_vol, track_particle, Kappa, dot_product); 

        checkTimeStep(geomParams, simParams, pos, u, c,
                      neighbours, nb_neighbours, acc_vol);

        double dt = simParams.dt;

        #pragma omp parallel for   
        for (int i = 0; i < simParams.nb_tot_part; i++){

            rho[i] += dt * drhodt[i];

            if (rho[i] < 0){
                cerr << "Rho negative (" << rho[i] << ") at t:" << simParams.t << endl;
                exit(EXIT_FAILURE);
            }

            for (int coord = 0; coord < 3; coord++){

                u[3*i + coord] += dt * dudt[3*i + coord];
                pos[3*i + coord] += dt * u[3*i + coord];
            }
        }
    }

    else if (scheme_integration == "RK22"){
    
        vector<double> u_half = u,
                       rho_half = rho,
                       pos_half = pos,
                       drhodt_half(nb_tot_part, 0.0),
                       dudt_half(3*nb_tot_part, 0.0),
                       acc_vol_half(3*nb_tot_part, 0.0);

        // First time step of RK22
       
        continuityEquation(simParams, neighbours, nb_neighbours, gradW, 
                           pos, u, drhodt, rho, mass);
        //printArray(drhodt,drhodt.size(), "drhodt");
        momentumEquation(geomParams, thermoParams, simParams, neighbours, nb_neighbours, gradW, 
                         W, viscosity, mass, dudt, rho, p, c, pos, u, type, colour, R, L, N, normal, acc_vol, track_particle, Kappa, dot_product); 
                        
        //checkTimeStep(geomParams, thermoParams, simParams, pos, u, c,
                     //neighbours, nb_neighbours);

        double theta = simParams.theta;
        double dt_half = simParams.dt/(2*theta);
        
        
                
        #pragma omp parallel for
        for (int n = 0; n < simParams.nb_tot_part; n++){

            rho_half[n] += dt_half * drhodt[n];
               
            if (rho_half[n] < 0){
                
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
                         W, viscosity, mass, dudt_half, rho_half, p, c, pos_half, u_half, type, colour, R, L, N, normal, acc_vol_half, track_particle, Kappa, dot_product); 


        checkTimeStep(geomParams, simParams, pos, u, c,
                      neighbours, nb_neighbours, acc_vol);

        double dt = simParams.dt;
        #pragma omp parallel for
        for (int n = 0; n < simParams.nb_tot_part; n++){
            
            rho[n] += dt * ((1-theta)*drhodt[n] + theta*drhodt_half[n]);
            
            for (int coord = 0; coord < 3; coord++){
                                
                pos[3 * n + coord] += dt * ((1-theta)*u[3 * n + coord] + theta*u_half[3 * n + coord]);
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

void checkTimeStep(GeomData geomParams,
                   SimulationData &simParams, 
                   vector<double> &pos,
                   vector<double> &u,
                   vector<double> &c,
                   vector<int> &neighbours,
                   vector<double> &nb_neighbours,
                   vector<double> &acc_vol){

    double alpha = simParams.alpha;
    double beta = simParams.beta;
    double h = geomParams.h;
    int nb_moving_part = simParams.nb_moving_part;
    int t = simParams.t;

    // Compute dt_f
    double dt_f = sqrt(h / simParams.acc_st_max);
    double dt_cv = 0;

    // Compute dt_cv
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
        for (int i = 0; i < nb_moving_part; i++){

            double c_a = c[i];
            max_b = 0; // always set to 0 in simulations
            double val = geomParams.h/(c_a + 0.6*(alpha*c_a + beta*max_b));
            min_a = (val < min_a) ? val : min_a;
        }

        dt_cv = min_a;
        double dt_final = min(0.4*dt_f, 0.25*dt_cv);
        string state_equation = simParams.state_equation;

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
