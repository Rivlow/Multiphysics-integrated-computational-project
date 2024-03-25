#include <stdio.h>
#include <vector>
#include "sorted_list.h"
using namespace std;

void gradW(size_t nb_moving_part,
           vector<vector<double>> &gradW_matrix,
           vector<double> &pos_arr,
           vector<vector<int>> &neighbours_matrix,
           double h, int Nx, int Ny, int Nz);

void setArtificialViscosity(int t, 
                            vector<vector<double>> &artificial_visc,
                            size_t nb_moving_part,
                            vector<double> &pos_arr,
                            vector<vector<int>> &neighbours_matrix,
                            vector<double> &rho_arr,
                            vector<double> &u_arr, 
                            double alpha, double beta, double gamma,
                            double c_0, double rho_0,
                            double R, double T,double M, 
                            double h,
                            string state_equation_chosen);

void continuityEquation(size_t nb_moving_part,
                        vector<double> &pos_arr,
                        vector<double> &u_arr,
                        vector<vector<int>> &neighbours_matrix,
                        vector<vector<double>> &gradW_matrix,
                        vector<double> &drhodt_arr,
                        vector<double> &rho_arr,
                        vector<double> &mass_arr, 
                        double h);

void momentumEquation(size_t nb_moving_part,
                      vector<vector<int>> &neighbours_matrix,
                      vector<double> &mass_arr,
                      vector<vector<double>> &gradW_matrix,
                      vector<double> &dudt_arr,
                      vector<vector<double>> &artificial_visc_matrix,
                      vector<double> &rho_arr,
                      vector<double> &p_arr, 
                      double rho_0, double c_0,
                      double gamma,
                      double R, double T, double M,
                      double g,
                      string state_equation_chosen);

double setSpeedOfSound(double rho, double rho_0, double c_0,
                       double gamma, 
                       string state_equation_chosen);

void setPressure(size_t nb_moving_part, 
                 vector<double> &p_arr,
                 vector<double> &rho_arr, 
                 double rho_0, double c_0, 
                 double R, double T, double M,
                 double gamma,
                 string state_equation_chosen);