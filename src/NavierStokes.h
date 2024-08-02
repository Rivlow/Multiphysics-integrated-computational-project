#ifndef GRADIENT_H
#define GRADIENT_H

#include <stdio.h>
#include <vector>
#include "structure.h"


using namespace std;

void computeGradW(GeomData &geomParams,    
                  SimulationData &simParams, 
                  vector<vector<double>> &gradW,
                  vector<vector<double>> &W,
                  vector<int> neighbours,
                  vector<double> nb_neighbours,
                  vector<double> pos);

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
                            vector<vector<double>> &viscosity,
                            vector<int> &neighbours,
                            vector<double> &nb_neighbours,
                            vector<double> &c,
                            vector<double> &pos,
                            vector<double> &rho,
                            vector<double> &u);

void continuityEquation(SimulationData &simParams, 
                        vector<int> &neighbours,
                        vector<double> &nb_neighbours,
                        vector<vector<double>> &gradW,
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
                      vector<vector<double>> &gradW,
                      vector<vector<double>> W,
                      vector<vector<double>> &viscosity,
                      vector<double> &mass,
                      vector<double> &dudt,
                      vector<double> &rho,
                      vector<double> &p, 
                      vector<double> &c,
                      vector<double> &pos,
                      vector<double> &u,
                      vector<double> type,
                      vector<double> &colour,
                      vector<double> &R,
                      vector<double> &L,
                      vector<double> &N,
                      vector<double> &normal,
                      vector<double> &F_vol,
                      vector<double> &track_particle,
                      vector<double> &Kappa,
                      vector<double> &dot_product);

#endif // GRADIENT_H
