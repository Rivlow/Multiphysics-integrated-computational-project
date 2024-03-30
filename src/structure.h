#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <stdio.h>
#include <vector>
#include <string>
#include <filesystem>

using namespace std;

struct SimulationData {

    int kappa;
    int nstepT;
    int nsave;
    double dt; 
    double s;
    double h; 
    vector<double> o;
    vector<double> L;
    vector<double> o_d;
    vector<double> L_d;
    vector<double> u_init;

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
};

#endif // STRUCTURE_H
