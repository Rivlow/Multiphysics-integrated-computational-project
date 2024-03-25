#include <stdio.h>
#include <vector>
#include <string>
using namespace std;

void initializeMass(vector<double> &rho_array,
                    vector<double> &mass_array,
                    double s);

void initializeRho(vector<double> &pos_array,
                   vector<double> &rho_array, 
                   size_t nb_moving_part,
                   double rho_moving,
                   double rho_fixed,
                   double rho_0, double c_0,
                   double M, double g, double R,
                   double T, double gamma,
                   std::string state_equation_chosen,
                   std::string state_initial_condition);

void initializeVelocity(vector<double> &u_array,
                        vector<double> &u_init,
                        size_t nb_moving_part);

void initializeViscosity(vector<vector<double>> &artificial_visc_matrix);