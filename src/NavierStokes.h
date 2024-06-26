#ifndef GRADIENT_H
#define GRADIENT_H

#include <stdio.h>
#include <vector>
#include "structure.h"


using namespace std;

void computeGradW(GeomData &geomParams,    
                   SimulationData &simParams, 
                   vector<double> &gradW,
                   vector<double> &W,
                   vector<int> neighbours,
                   vector<int> nb_neighbours,
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
                            vector<double> &viscosity,
                            vector<int> &neighbours,
                            vector<int> &nb_neighbours,
                            vector<double> &c,
                            vector<double> &pos,
                            vector<double> &rho,
                            vector<double> &u);

void continuityEquation(SimulationData &params,
                        vector<int> &neighbours,
                        vector<int> &nb_neighbours,
                        vector<double> &gradW,
                        vector<double> &pos,
                        vector<double> &u,
                        vector<double> &drhodt,
                        vector<double> &rho,
                        vector<double> &mass);

void momentumEquation(GeomData &geomParams,    
                      ThermoData &thermoParams,
                      SimulationData &simParams, 
                      vector<int> &neighbours,
                      vector<int> &nb_neighbours,
                      vector<double> &gradW,
                      vector<double> W,
                      vector<double> &viscosity,
                      vector<double> &mass,
                      vector<double> &dudt,
                      vector<double> &rho,
                      vector<double> &p, 
                      vector<double> &c,
                      vector<double> &pos,
                      vector<double> &u,
                      vector<int> type,
                      vector<int> &track_particle);

#endif // GRADIENT_H
