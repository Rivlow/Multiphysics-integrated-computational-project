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
#include "integration.h"
#include "structure.h"
#include "gradient.h"
#include "tools.h"
#include <omp.h>

void Euler( SimulationData& params,
            int t,
            vector<double> &pos,
            vector<double> &u,
            vector<double> &rho,
            vector<double> &drhodt,
            vector<double> &c,
            vector<double> &p,
            vector<double> &dudt,
            vector<double> &mass,
            vector<vector<double>> &artificial_visc_matrix,
            vector<vector<double>> &gradW_matrix,
            vector<vector<int>> &neighbours_matrix){


    string schemeIntegration = params.schemeIntegration;
    double theta = params.theta;
    double dt = params.dt;

    if (schemeIntegration == "RK22"){
        dt = params.dt/(2*theta);
    }

    if (schemeIntegration == "Euler"){
        dt = params.dt;
    }

    // Compute D(rho)/Dt for all particles
    continuityEquation(params, neighbours_matrix, gradW_matrix, 
                    pos, u, drhodt, rho, mass); 

    // Compute D(u)/Dt for all particles
    momentumEquation(params, t, neighbours_matrix, gradW_matrix, artificial_visc_matrix,
                    mass, dudt, rho, p, c, pos, u); 

    
    int size_pos = pos.size()/3;
    #pragma omp parallel for   
    for (int n = 0; n < size_pos; n++){

        rho[n] += dt * drhodt[n];

        for (int cord = 0; cord < 3; cord++){

            pos[3 * n + cord] += dt * u[3 * n + cord];
            u[3 * n + cord] += dt * dudt[3 * n + cord];
        }
    }
}

void RK22( SimulationData& params,
          int t,
          vector<double> &pos,
          vector<double> &u,
          vector<double> &rho,
          vector<double> &drhodt,
          vector<double> &c,           
          vector<double> &p,
          vector<double> &dudt,
          vector<double> &mass,
          vector<vector<double>> &artificial_visc_matrix,
          vector<vector<double>> &gradW_matrix,
          vector<vector<int>> &neighbours_matrix){

    double dt = params.dt;
    double theta = params.theta;  
    string schemeIntegration = params.schemeIntegration;
    int nb_tot_part = pos.size()/3;

    vector<double>  u_half = u,
                    rho_half = rho,
                    pos_half = pos,
                    drhodt_half(nb_tot_part,0.0),
                    dudt_half(3*nb_tot_part,0.0);

    Euler(params, t, pos_half, u_half, rho_half, drhodt_half, c, p, dudt_half, 
              mass, artificial_visc_matrix, gradW_matrix, neighbours_matrix);

    // Compute D(rho)/Dt for all particles
    continuityEquation(params, neighbours_matrix, gradW_matrix, 
                    pos_half, u_half, drhodt_half, rho_half, mass); 

    // Compute D(u)/Dt for all particles
    momentumEquation(params, t, neighbours_matrix, gradW_matrix, artificial_visc_matrix,
                    mass, dudt_half, rho_half, p, c, pos_half, u_half); 

    int size_pos = pos.size()/3;

    #pragma omp parallel for
    for (int n = 0; n < size_pos; n++){

        rho[n] += dt * ((1-theta)*drhodt[n] + theta*drhodt_half[n]);

        for (int cord = 0; cord < 3; cord++){

            double u_temp = u[3 * n + cord] + dt*dudt_half[3 * n + cord];
            pos[3 * n + cord] += dt * ((1-theta)*u[3 * n + cord] + theta*u_temp);
            u[3 * n + cord] += dt *  ((1-theta)*dudt[3 * n + cord] + theta*dudt_half[3 * n + cord]);

        }
    } 
}

void updateVariables( SimulationData& params,
                     int t,
                     vector<double> &pos,
                     vector<double> &u,
                     vector<double> &rho,
                     vector<double> &drhodt,
                     vector<double> &c,
                     vector<double> &p,
                     vector<double> &dudt,
                     vector<double> &mass,
                     vector<vector<double>> &artificial_visc_matrix,
                     vector<vector<double>> &gradW_matrix,
                     vector<vector<int>> &neighbours_matrix){

    bool PRINT = params.PRINT;
    string schemeIntegration = params.schemeIntegration;
    int nb_tot_part = pos.size()/3;

    if (schemeIntegration == "Euler"){
    
        Euler(params, t, pos, u, rho, drhodt, c, p, dudt, mass, artificial_visc_matrix, gradW_matrix, neighbours_matrix);
    }

    if (schemeIntegration == "RK22"){

        vector<double>  u_half = u,
                        rho_half = rho,
                        pos_half = pos,
                        drhodt_half(nb_tot_part,0.0),
                        dudt_half(3*nb_tot_part,0.0);

        RK22(params, t, pos, u, rho, drhodt, c, p, dudt, mass, artificial_visc_matrix, gradW_matrix, neighbours_matrix);

    }

    if (PRINT){
        cout << "update passed" << endl;
    }
}
