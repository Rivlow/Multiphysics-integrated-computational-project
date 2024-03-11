#include <stdio.h>
#include <vector>
#include <string.h>
using namespace std;


void initializeMass(vector<double> &rho_arr, double &s, vector<double> &mass_arr);

void initializeRho(vector<double> &rho_arr,double &rho, double &rho_0, double &c_0, 
                   double &M, double &g, double &H, double &R, double &T, double &gamma, string &state_equation_chosen, 
                   string &state_initial_condition);

void initializeVelocity(vector<double> &u_arr, vector<double> &u_init);