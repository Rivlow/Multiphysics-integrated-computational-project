#include <stdio.h>
#include <vector>
#include <string>
#include <limits>
#include <iostream>

#include "initialize.h"
#include "structure.h"
#include "NavierStokes.h"
#include "tools.h"

using namespace std;


void initMass(GeomData &geomParams,
              SimulationData &simParams, 
              vector<double> &rho,
              vector<double> &mass){

    double s = geomParams.s;  
    bool PRINT = simParams.PRINT; 
    int dim = simParams.dimension;
    double V = 0.0;

    if (dim == 3)
        V = s * s * s;
    else
        V = s*s;

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
    double g = (simParams.is_gravity)? -9.81 : 0;

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
                rho[i] = (i < nb_moving_part) ? rho_0 * pow(1 + rho_0 
                                * g * pos[3 * i + 2] / B, 1/gamma) : rho_fixed;  
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

    if (PRINT)
       cout << "initVelocity passed" << endl;
    
}



void initKernelCoef(GeomData &geomParams, 
                    SimulationData &simParams){

    double h = geomParams.h;

    simParams.cubic_kernel_coef = (simParams.dimension == 2)? 15.0 /( 7.0 * M_PI * h * h ) : 3.0 / (2.0 * M_PI * h * h * h);
    simParams.adh_kernel_coef = (simParams.dimension == 2)? 16/(4* M_PI*pow(h, 2.25)): 0.0007/pow(h,3.25);;
    simParams.coh_kernel_coef = (simParams.dimension == 2) ? 40/(M_PI*h*h*h*h*h*h*h*h) : 32/(M_PI*h*h*h*h*h*h*h*h*h);
    simParams.quintic_kernel_coef = (simParams.dimension == 2)? 7.0 /( 4.0 * M_PI * h * h ) : 0;


}
