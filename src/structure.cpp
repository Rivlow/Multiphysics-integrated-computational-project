#include <stdio.h>
#include <string.h>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include "nlohmann/json.hpp"
#include "structure.h"
#include "generate_particle.h"

using json = nlohmann::json;
using namespace std;

//extern GeomData geomParams;
//extern ThermoData thermoParams;
//extern SimulationData simParams;


void initializeStruct(json data,
                      string state_equation, 
                      string state_initial_condition, 
                      string schemeIntegration,
                      string kernel,
                      GeomData &geomParams,
                      SimulationData &simParams,
                      ThermoData &thermoParams){

    
    geomParams.kappa = data["simulation"]["kappa"];
    geomParams.s = data["simulation"]["s"];
    geomParams.h = 1.2 * data["simulation"]["s"].get<double>();
    geomParams.o_d = data["domain"]["o_d"].get<vector<double>>();
    geomParams.L_d = data["domain"]["L_d"].get<vector<double>>();
    geomParams.matrix_long = data["domain"]["matrix_long"].get<vector<vector<double>>>();
    geomParams.matrix_orig = data["domain"]["matrix_orig"].get<vector<vector<double>>>();
    geomParams.vector_type = data["domain"]["vector_type"].get<vector<int>>();
    geomParams.sphere_do = data["domain"]["sphere"]["do"].get<vector<int>>();
    geomParams.radius = data["domain"]["sphere"]["radius"].get<vector<double>>();
    geomParams.xyz_init = data["post_process"]["xyz_init"].get<vector<double>>();
    geomParams.xyz_end = data["post_process"]["xyz_end"].get<vector<double>>();
    geomParams.post_process_do = data["post_process"]["do"];
    geomParams.Nx = int(geomParams.L_d[0] / (geomParams.kappa * geomParams.h));
    geomParams.Ny = int(geomParams.L_d[1] / (geomParams.kappa * geomParams.h));
    geomParams.Nz = int(geomParams.L_d[2] / (geomParams.kappa * geomParams.h));
    geomParams.following_part_bool = data["following_part"]["part"];
    geomParams.following_part_min = data["following_part"]["min"];
    geomParams.following_part_max = data["following_part"]["max"];
    geomParams.following_part_part = data["following_part"]["particle"];
    geomParams.following_part_p = data["following_part"]["pressure"];
    geomParams.following_part_rho = data["following_part"]["rho"];
    geomParams.following_part_pos = data["following_part"]["position"].get<vector<bool>>();
    geomParams.following_part_u = data["following_part"]["velocity"].get<vector<bool>>();
    cout << "geomParams initialized" << endl;

    thermoParams.c_0 = data["thermo"]["c_0"];
    thermoParams.rho_moving = data["thermo"]["rho_moving"];
    thermoParams.rho_fixed = data["thermo"]["rho_fixed"];
    thermoParams.rho_0 = data["thermo"]["rho_0"];
    thermoParams.M = data["thermo"]["M"];
    thermoParams.T = data["thermo"]["T"];
    thermoParams.gamma = data["thermo"]["gamma"];
    thermoParams.R = data["thermo"]["R"];
    thermoParams.sigma = data["thermo"]["sigma"];
    cout << "thermoParams initialized" << endl;


    simParams.dimension = data["simulation"]["dimension"];
    simParams.nstepT = data["simulation"]["nstepT"];
    simParams.nsave = data["simulation"]["nsave"];
    simParams.dt = data["simulation"]["dt"];
    simParams.theta = data["simulation"]["theta"];
    simParams.alpha = data["simulation"]["alpha"];
    simParams.beta = data["simulation"]["beta"];
    simParams.alpha_st = data["simulation"]["alpha_st"];
    simParams.beta_adh = data["simulation"]["beta_adh"];
    simParams.schemeIntegration = schemeIntegration;
    simParams.u_init = data["thermo"]["u_init"].get<vector<double>>();
    simParams.state_equation = state_equation;
    simParams.state_initial_condition = state_initial_condition;
    simParams.kernel = kernel;
    simParams.is_gravity = data["forces"]["gravity"];
    simParams.is_surface_tension_1 = data["forces"]["surface_tension_1"];
    simParams.is_surface_tension_2 = data["forces"]["surface_tension_2"];
    simParams.is_adhesion = data["forces"]["adhesion"];
    simParams.PRINT = data["condition"]["print_debug"];
    simParams.comparison_algorithm  = data["simulation"]["comparison_algorithm"];
    simParams.nb_moving_part = evaluateNumberParticles(geomParams),
    simParams.nb_tot_part = 0;
    simParams.t = 0;
    simParams.acc_st_max = 0.0;
    simParams.cubic_kernel_coef = 0.0;
    simParams.adh_kernel_coef = 0.0;
    simParams.coh_kernel_coef = 0.0;
    simParams.quintic_kernel_coef = 0.0;
    cout << "simParams initialized" << endl;

}
