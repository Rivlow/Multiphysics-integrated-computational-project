#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void Euler(const SimulationData& params,
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
            vector<vector<int>> &neighbours_matrix);

void RK22(const SimulationData& params,
           vector<double> &rho_array,
           vector<double> &pos_array,
           vector<double> &u_array,
           vector<double> drhodt_array,
           vector<double> drhodt_array_half,
           vector<double> dudt_array,
           vector<double> dudt_array_half);

void updateVariables(const SimulationData& params,
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
            vector<vector<int>> &neighbours_matrix);
