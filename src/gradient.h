#ifndef GRADIENT_H
#define GRADIENT_H

#include <stdio.h>
#include <vector>
#include "structure.h"


using namespace std;

void gradW(GeomData &geomParams,    
           SimulationData &simParams, 
           vector<vector<double>> &gradW_matrix,
           vector<vector<double>> &W_matrix,
           vector<int> &neighbours,
           vector<double> &nb_neighbours,
           vector<double> &pos);

void setSpeedOfSound(GeomData &geomParams,    
                     ThermoData &thermoParams,
                     SimulationData &simParams, 
                     vector<double> &c,
                     vector<double> &rho);

void setPressure(GeomData &geomParams,    
                 ThermoData &thermoParams,
                 SimulationData &simParams, 
                 vector<double> &p,
                 vector<double> &rho);

void setArtificialViscosity(GeomData &geomParams,    
                            ThermoData &thermoParams,
                            SimulationData &simParams, 
                            vector<vector<double>> &pi_matrix,
                            vector<int> &neighbours,
                            vector<double> &nb_neighbours,
                            vector<double> &c,
                            vector<double> &pos,
                            vector<double> &rho,
                            vector<double> &u);

void continuityEquation(SimulationData &params,
                        vector<int> &neighbours,
                        vector<double> &nb_neighbours,
                        vector<vector<double>> &gradW_matrix,
                        vector<double> &pos,
                        vector<double> &u,
                        vector<double> &drhodt,
                        vector<double> &rho,
                        vector<double> &mass);

void momentumEquation(GeomData &geomParams,    
                      ThermoData &thermoParams,
                      SimulationData &simParams, 
                      vector<int> &neighbours,
                      vector<double> &nb_neighbours,
                      vector<int> &track_surface,
                      vector<double> &N_smoothed,
                      vector<vector<double>> &gradW_matrix,
                      vector<vector<double>> W_matrix,
                      vector<vector<double>> &pi_matrix,
                      vector<double> &mass,
                      vector<double> &dudt,
                      vector<double> &rho,
                      vector<double> &p, 
                      vector<double> &c,
                      vector<double> &pos,
                      vector<double> &u,
                      vector<double> type);

#endif // GRADIENT_H