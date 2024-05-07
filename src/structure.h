#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <stdio.h>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include "nlohmann/json.hpp"

using namespace std;

struct GeomData {
    
    int kappa;
    double s;
    double h; 
    vector<double> o_d;
    vector<double> L_d;
    vector<double> xyz_init;
    vector<double> xyz_end;
    bool post_process_do;
    int Nx;
    int Ny;
    int Nz;
    vector<vector<double>> matrix_long;
    vector<vector<double>> matrix_orig;
    vector<int> vector_type; 
};

struct ThermoData {
 
    double alpha;
    double beta;
    double c_0;
    double rho_moving;
    double rho_fixed;
    double rho_0;
    double M;
    double T;
    double gamma;
    double R ; // [J/(K.mol)]
    double g; // [m/s²]

};

struct SimulationData {

    int nstepT;
    int nsave;
    double dt; 
    double theta;
    string schemeIntegration;
    vector<double> u_init;
    string state_equation; 
    string state_initial_condition;
    bool PRINT;
    int nb_moving_part;
    int nb_part;
    int t;
 
};

#endif // STRUCTURE_H
