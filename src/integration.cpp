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

void Euler(GeomData &geomParams,    
           ThermoData &thermoParams,
           SimulationData &simParams, 
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
            vector<vector<int>> &neighbours_matrix,
            vector<double> &nb_neighbours){


    string schemeIntegration = simParams.schemeIntegration;
    double theta = simParams.theta;
    double dt = simParams.dt;

    if (schemeIntegration == "RK22"){
        dt = simParams.dt/(2*theta);
    }

    if (schemeIntegration == "Euler"){
        dt = simParams.dt;
    }

    // Compute D(rho)/Dt for all particles
    continuityEquation(simParams, neighbours_matrix, nb_neighbours, gradW_matrix, 
                    pos, u, drhodt, rho, mass); 

    // Compute D(u)/Dt for all particles
    momentumEquation(geomParams, thermoParams, simParams, t, neighbours_matrix, nb_neighbours, gradW_matrix, artificial_visc_matrix,
                    mass, dudt, rho, p, c, pos, u); 

    
    int nb_fixed_part = simParams.nb_fixed_part;

    #pragma omp parallel for   
    for (int n = 0; n < nb_fixed_part; n++){

        rho[n] += dt * drhodt[n];

        for (int cord = 0; cord < 3; cord++){

            pos[3 * n + cord] += dt * u[3 * n + cord];
            u[3 * n + cord] += dt * dudt[3 * n + cord];
        }
    }
}

void RK22(GeomData &geomParams,    
          ThermoData &thermoParams,
          SimulationData &simParams, 
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
          vector<vector<int>> &neighbours_matrix,
          vector<double> &nb_neighbours){

    double dt = simParams.dt;
    double theta = simParams.theta;  
    string schemeIntegration = simParams.schemeIntegration;
    int nb_tot_part = pos.size()/3;

    vector<double>  u_half = u,
                    rho_half = rho,
                    pos_half = pos,
                    drhodt_half(nb_tot_part,0.0),
                    dudt_half(3*nb_tot_part,0.0);

    Euler(geomParams, thermoParams, simParams, t, pos_half, u_half, rho_half, drhodt_half, c, p, dudt_half, 
              mass, artificial_visc_matrix, gradW_matrix, neighbours_matrix, nb_neighbours);

    // Compute D(rho)/Dt for all particles
    continuityEquation(simParams, neighbours_matrix, nb_neighbours, gradW_matrix, 
                    pos_half, u_half, drhodt_half, rho_half, mass); 

    // Compute D(u)/Dt for all particles
    momentumEquation(geomParams, thermoParams, simParams, t, neighbours_matrix, nb_neighbours, gradW_matrix, artificial_visc_matrix,
                    mass, dudt_half, rho_half, p, c, pos_half, u_half); 

    int nb_fixed_part = simParams.nb_fixed_part;

    #pragma omp parallel for
    for (int n = 0; n < nb_fixed_part; n++){

        rho[n] += dt * ((1-theta)*drhodt[n] + theta*drhodt_half[n]);

        for (int cord = 0; cord < 3; cord++){

            double u_temp = u[3 * n + cord] + dt*dudt_half[3 * n + cord];
            pos[3 * n + cord] += dt * ((1-theta)*u[3 * n + cord] + theta*u_temp);
            u[3 * n + cord] += dt *  ((1-theta)*dudt[3 * n + cord] + theta*dudt_half[3 * n + cord]);

        }
    } 
}

void updateVariables(GeomData &geomParams,    
                     ThermoData &thermoParams,
                     SimulationData &simParams, 
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
                     vector<vector<int>> &neighbours_matrix,
                     vector<double> &nb_neighbours){

    bool PRINT = simParams.PRINT;
    string schemeIntegration = simParams.schemeIntegration;
    int nb_tot_part = pos.size()/3;

    if (schemeIntegration == "Euler"){
    
        Euler(geomParams, thermoParams, simParams, t, pos, u, rho, drhodt, c, p, dudt, mass, artificial_visc_matrix, gradW_matrix, neighbours_matrix, nb_neighbours);
    }

    if (schemeIntegration == "RK22"){

        vector<double>  u_half = u,
                        rho_half = rho,
                        pos_half = pos,
                        drhodt_half(nb_tot_part,0.0),
                        dudt_half(3*nb_tot_part,0.0);

        RK22(geomParams, thermoParams, simParams, t, pos, u, rho, drhodt, c, p, dudt, mass, artificial_visc_matrix, gradW_matrix, neighbours_matrix, nb_neighbours);

    }

    if (PRINT){
        cout << "update passed" << endl;
    }
}
