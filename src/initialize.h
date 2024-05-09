#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <stdio.h>
#include <vector>
#include <string>

#include "tools.h"
#include "structure.h"


using namespace std;

void initMass(GeomData &geomParams,
                    SimulationData &simParams, 
                    vector<double> &rho,
                    vector<double> &mass);

void initRho(ThermoData &thermoParams,
                   SimulationData &simParams,
                   vector<double> &pos,
                   vector<double> &rho);

void initVelocity(ThermoData &thermoParams,
                        SimulationData &simParams, 
                        vector<double> &u);

void initViscosity(SimulationData &simParams, 
                         vector<vector<double>> &pi_matrix);


#endif // INITIALIZE_H