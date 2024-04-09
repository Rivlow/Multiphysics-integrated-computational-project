#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <stdio.h>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include "nlohmann/json.hpp"

using namespace std;


struct DomainParams {
    string shape;
    vector<string> walls_used;
    double particle_layers;

};

struct SimulationData {
    int kappa;
    int nstepT;
    int nsave;
    double dt; 
    double theta;
    string schemeIntegration;
    double s;
    double h; 
    vector<double> o;
    vector<double> L;
    vector<double> o_d;
    vector<double> L_d;
    vector<double> u_init;
    int Nx;
    int Ny;
    int Nz;

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

    string state_equation; 
    string state_initial_condition;
    bool PRINT;

    int nb_moving_part;

    DomainParams domainParams;

    bool walls_used(string wall)  {
        return find(domainParams.walls_used.begin(), domainParams.walls_used.end(), wall) != domainParams.walls_used.end();
    }
};

#endif // STRUCTURE_H
