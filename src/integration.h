#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void Euler(const double dt, 
           const size_t nb_tot_part,
           vector<double> &rho_array,
           vector<double> &pos_array,
           vector<double> &u_array,
           vector<double> drhodt_array,
           vector<double> dudt_array);

void runchKutta(const double dt, 
           const size_t nb_tot_part,
           const double theta,
           vector<double> &rho_array,
           vector<double> &pos_array,
           vector<double> &u_array,
           vector<double> drhodt_array,
           vector<double> drhodt_array_half,
           vector<double> dudt_array,
           vector<double> dudt_array_half);