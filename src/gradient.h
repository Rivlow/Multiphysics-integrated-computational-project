#include <stdio.h>
#include <vector>
#include "sorted_list.h"
using namespace std;

void gradW(vector<vector<double>> &gradW_matrix,
           vector<vector<int>> &neighbours_matrix,
           vector<double> &pos_array,
           size_t nb_moving_part,
           double h, int Nx, int Ny, int Nz, 
           const bool PRINT);

void setArtificialViscosity(vector<vector<double>> &artificial_visc_matrix,
                            vector<vector<int>> &neighbours_matrix,
                            vector<double> &pos_array,
                            vector<double> &rho_array,
                            vector<double> &u_array, 
                            size_t nb_moving_part,
                            int t, 
                            double alpha, double beta, double gamma,
                            double c_0, double rho_0,
                            double R, double T,double M, 
                            double h,
                            string state_equation_chosen, 
                            const bool PRINT);

void continuityEquation(vector<vector<int>> &neighbours_matrix,
                        vector<vector<double>> &gradW_matrix,
                        vector<double> &pos_array,
                        vector<double> &u_array,
                        vector<double> &drhodt_array,
                        vector<double> &rho_array,
                        vector<double> &mass_array, 
                        size_t nb_moving_part,
                        double h, 
                        const bool PRINT);

void momentumEquation(vector<vector<int>> &neighbours_matrix,
                      vector<vector<double>> &gradW_matrix,
                      vector<vector<double>> &artificial_visc_matrix,
                      vector<double> &mass_array,
                      vector<double> &dudt_array,
                      vector<double> &rho_array,
                      vector<double> &p_array, 
                      size_t nb_moving_part,
                      double rho_0, double c_0,
                      double gamma,
                      double R, double T, double M,
                      double g,
                      string state_equation_chosen, 
                      const bool PRINT);

double setSpeedOfSound(double rho, double rho_0, double c_0,
                       double gamma, 
                       string state_equation_chosen);

void setPressure(vector<double> &p_array,
                 vector<double> &rho_array, 
                 size_t nb_moving_part,
                 double rho_0, double c_0, 
                 double R, double T, double M,
                 double gamma,
                 string state_equation_chosen, 
                 const bool PRINT);