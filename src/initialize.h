#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <stdio.h>
#include <vector>
#include <string>

#include "tools.h"
#include "structure.h"


using namespace std;

void initializeMass(GeomData &geomParams,
                    SimulationData &simParams, 
                    vector<double> &rho,
                    vector<double> &mass);

void initializeRho(ThermoData &thermoParams,
                   SimulationData &simParams,
                   vector<double> &pos,
                   vector<double> &rho);

void initializeVelocity(ThermoData &thermoParams,
                        SimulationData &simParams, 
                        vector<double> &u);

void initializeViscosity(SimulationData &simParams, 
                         vector<vector<double>> &pi_matrix);

void checkTimeStep(GeomData &geomParams,    
                   ThermoData &thermoParams,
                   SimulationData &simParams, 
                   int t,
                   vector<double> pos,
                   vector<double> c,
                   vector<vector<int>> &neighbours_matrix,
                   vector<double> &nb_neighbours,
                   vector<vector<double>> &pi_matrix);

#endif // INITIALIZE_H