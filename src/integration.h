#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void Euler( SimulationData& params,
            int t,
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
            vector<double> &nb_neighbours);

void RK22( SimulationData& params,
          int t,
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
          vector<double> &nb_neighbours);

void updateVariables( SimulationData& params,
            int t,
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
            vector<double> &nb_neighbours);
