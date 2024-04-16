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
    vector<double> o;
    vector<double> L;
    vector<double> o_d;
    vector<double> L_d;
    int Nx;
    int Ny;
    int Nz;
    
    vector<vector<double>> matrixLong;
    vector<vector<double>> matrixOrig;
    vector<int> vectorType; 

    
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
    double R = 8.314; // [J/(K.mol)]
    double g = -9.81; // [m/s²]

};

struct SimulationData {

    int nstepT;
    int nsave;
    double dt; 
    double theta;
    string schemeIntegration;
    vector<string> data_store;
    int data_init;
    int data_end;
    bool data_do;
    
    vector<double> u_init;

    string state_equation; 
    string state_initial_condition;
    bool PRINT;

    int nb_moving_part;

 
};

#endif // STRUCTURE_H
