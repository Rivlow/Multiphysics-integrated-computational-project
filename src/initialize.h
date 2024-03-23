#include <stdio.h>
#include <vector>
#include <string>
using namespace std;


void initializeMass(vector<double> &rho_arr, double &s, vector<double> &mass_arr);

void initializeRho(int &nb_moving_part, vector<double> &pos_arr, vector<double> &rho_arr, double &rho_moving, double &rho_fixed, double &rho_0, double &c_0, 
                   double &M, double &g, double &R, double &T, double &gamma, std::string &state_equation_chosen, 
                   std::string &state_initial_condition);

void initializeVelocity(int &nb_moving_part, vector<double> &u_arr, vector<double> &u_init);

void initializeViscosity(vector<vector<double>> &artificial_visc_matrix);