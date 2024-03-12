#include <stdio.h>
#include <vector>
#include <string>
using namespace std;


void initializeMass(vector<double> &rho_arr, double &s, vector<double> &mass_arr);

void initializeRho(vector<double> &pos_arr, vector<double> &rho_arr,double &rho, double &rho_0, double &c_0, 
                   double &M, double &g, double &R, double &T, double &gamma, std::string &state_equation_chosen, 
                   std::string &state_initial_condition);

void initializeVelocity(vector<double> &u_arr, vector<double> &u_init);