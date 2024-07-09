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
                     vector<double> &viscosity,
                     vector<double> &gradW,
                     vector<double> &W,
                     vector<int> &neighbours,
                     vector<double> &nb_neighbours,
                     vector<double> type,
                     vector<double> &track_particle);

void checkTimeStep(GeomData &geomParams,    
                   ThermoData &thermoParams,
                   SimulationData &simParams, 
                   vector<double> &pos,
                   vector<double> &u,
                   vector<double> &c,
                   vector<int> &neighbours,
                   vector<double> &nb_neighbours);
