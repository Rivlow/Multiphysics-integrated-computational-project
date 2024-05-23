#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void Euler(GeomData &geomParams,    
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
            vector<vector<double>> &artificial_visc_matrix,
            vector<vector<double>> &gradW_matrix,
            vector<vector<int>> &neighbours_matrix,
            vector<double> &nb_neighbours,
            vector<double> &normal_print);

void RK22(GeomData &geomParams,    
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
          vector<vector<double>> &artificial_visc_matrix,
          vector<vector<double>> &gradW_matrix,
          vector<vector<int>> &neighbours_matrix,
          vector<double> &nb_neighbours,
          vector<double> &normal_print);

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
                     vector<vector<double>> &artificial_visc_matrix,
                     vector<vector<double>> &gradW_matrix,
                     vector<vector<int>> &neighbours_matrix,
                     vector<double> &nb_neighbours,
                     vector<double> &normal_print);
