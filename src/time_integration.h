#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void updateVariables(GeomData &geomParams,    
                     ThermoData &thermoParams,
                     SimulationData &simParams, 
                     vector<double> &pos,
                     vector<double> &u,
                     vector<double> &rho,
                     vector<double> &drhodt,
                     vector<double> &c,
                     vector<double> &p,
                     vector<double> &dudt,
                     vector<double> &mass,
                     vector<vector<double>> &viscosity,
                     vector<vector<double>> &gradW,
                     vector<vector<double>> &W,
                     vector<int> &neighbours,
                     vector<double> &nb_neighbours,
                     vector<double> type,
                     vector<double> &colour,
                     vector<double> &R,
                     vector<double> &L,
                     vector<double> &N,
                     vector<double> &normal,
                     vector<double> &acc_vol,
                     vector<double> &track_particle,
                     vector<double> &Kappa,
                     vector<double> &dot_product);

void checkTimeStep(GeomData geomParams,
                   SimulationData &simParams, 
                   vector<double> &pos,
                   vector<double> &u,
                   vector<double> &c,
                   vector<int> &neighbours,
                   vector<double> &nb_neighbours,
                   vector<double> &acc_vol);
