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
#include <omp.h>

void Euler(const SimulationData& params,
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

    // Compute D(rho)/Dt for all particles
    continuityEquation(params, neighbours_matrix, gradW_matrix, 
                    pos, u, drhodt, rho, mass); 

    // Compute D(u)/Dt for all particles
    momentumEquation(params, t, neighbours_matrix, gradW_matrix, artificial_visc_matrix,
                    mass, dudt, rho, p, c, pos, u); 

    
    #pragma omp parallel for   
    for (size_t idx = 0; idx < pos.size()/3; idx++){

            rho[idx] += dt * drhodt[idx];

            for (size_t cord = 0; cord < 3; cord++)
            {
                pos[3 * idx + cord] += dt * u[3 * idx + cord];
                u[3 * idx + cord] += dt * dudt[3 * idx + cord];
            }
        }
}

void RK22(const SimulationData& params,
          vector<double> &rho_array,
          vector<double> &pos_array,
          vector<double> &u_array,
          vector<double> drhodt_array,
          vector<double> drhodt_array_half,
          vector<double> dudt_array,
          vector<double> dudt_array_half){

    double dt = params.dt;
    double theta = params.theta;  

    #pragma omp parallel for
    for (size_t pos = 0; pos < pos_array.size()/3; pos++){

        rho_array[pos] += dt * ((1-theta)*drhodt_array[pos] + theta*drhodt_array_half[pos]);

        for (size_t cord = 0; cord < 3; cord++)
        {
            double u_temp = u_array[3 * pos + cord] + dt*dudt_array_half[3 * pos + cord];
            pos_array[3 * pos + cord] += dt * ((1-theta)*u_array[3 * pos + cord] + theta*u_temp);
            u_array[3 * pos + cord] += dt *  ((1-theta)*dudt_array[3 * pos + cord] + theta*dudt_array_half[3 * pos + cord]);

        }
    } 
}

void updateVariables(const SimulationData& params,
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
    bool PRINT = params.PRINT;
    double theta = params.theta;  
    string schemeIntegration = params.schemeIntegration;
    int nb_tot_part = pos.size()/3;

    if (schemeIntegration == "Euler"){
    
        Euler(params, t, pos, u, rho, drhodt, c, p, dudt, mass, artificial_visc_matrix, gradW_matrix, neighbours_matrix);
    }

    if (schemeIntegration == "RK22"){

        vector<double>  u_array_half = u,
                        rho_array_half = rho,
                        pos_array_half = pos,
                        drhodt_half(nb_tot_part,0.0),
                        dudt_half(3*nb_tot_part,0.0);

        Euler(params, t, pos, u, rho, drhodt, c, p, dudt, mass, artificial_visc_matrix, gradW_matrix, neighbours_matrix);

        // Compute D(rho)/Dt for all particles
        continuityEquation(params, neighbours_matrix, gradW_matrix, 
                        pos, u, drhodt, rho, mass); 

        // Compute D(u)/Dt for all particles
        momentumEquation(params, t, neighbours_matrix, gradW_matrix, artificial_visc_matrix,
                        mass, dudt, rho, p, c, pos, u); 

        //printArray(dudt_half,dudt_half.size(),"dudt_half");
        RK22(params, rho, pos, u, drhodt, drhodt_half, dudt, dudt_half);
    }

    if (PRINT){
        cout << "update passed" << endl;
    }
}
