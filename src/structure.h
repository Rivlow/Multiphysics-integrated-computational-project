#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <stdio.h>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include "nlohmann/json.hpp"

using json = nlohmann::json; 
using namespace std;

struct GeomData {
    
    int kappa;
    double s;
    double h; 
    vector<double> o_d;
    vector<double> L_d;
    vector<vector<double>> matrix_long;
    vector<vector<double>> matrix_orig;
    vector<int> vector_type; 
    vector<int> sphere_do;
    vector<double> radius;
    vector<double> xyz_init;
    vector<double> xyz_end;
    bool post_process_do;
    int Nx;
    int Ny;
    int Nz;
    bool following_part_bool;
    bool following_part_min;
    bool following_part_max;
    int following_part_part;
    bool following_part_p;
    bool following_part_rho;
    vector<bool> following_part_pos;
    vector<bool> following_part_u;

};



struct ThermoData {
    double c_0;
    double rho_moving;
    double rho_fixed;
    double rho_0;
    double M;
    double T;
    double gamma;
    double R ; // [J/(K.mol)]
    double sigma;
};

struct SimulationData {

    int dimension;
    int nstepT;
    int nsave;
    double dt; 
    double theta;
    double alpha;
    double beta;
    double alpha_st;
    double beta_adh;
    string schemeIntegration;
    vector<double> u_init;
    string state_equation; 
    string state_initial_condition;
    string kernel;
    bool is_gravity;
    bool is_surface_tension_1;
    bool is_surface_tension_2;
    bool is_adhesion;
    bool PRINT;
    bool comparison_algorithm;
    int nb_moving_part;
    int nb_tot_part;
    int t;
    double acc_st_max;
    double cubic_kernel_coef;
    double adh_kernel_coef;
    double coh_kernel_coef;
    double quintic_kernel_coef;

};

void initializeStruct(json data,
                      string state_equation, 
                      string state_initial_condition, 
                      string schemeIntegration,
                      string kernel,
                      GeomData &geomParams,
                      SimulationData &simParams,
                      ThermoData &thermoParams);


#endif // STRUCTURE_H
