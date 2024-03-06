#include <stdio.h>
#include <vector>
#include "sorted_list.h"
using namespace std;

void gradW( vector<double> &part_pos,  vector<vector<unsigned>> &neighbours_matrix, vector<vector<double>> &gradW_matrix, 
           double L[3],  double &h,  unsigned &Nx,  unsigned &Ny,  unsigned &Nz);

double setArtificialViscosity(vector<vector<double>> &artificial_visc,  vector<double> &part_pos,  vector<vector<unsigned>> &neighbours_matrix, vector<double> &u_arr, 
                             double &c_ab,  double &rho_ab,  double &alpha,  double &beta,  double &h);

void continuityEquation( vector<double> &part_pos,  vector<double> &u_arr,  vector<vector<unsigned>> &neighbours_matrix, 
                         vector<vector<double>> &gradW_matrix, vector<double> &drhodt_arr, vector<double> &rho_arr, vector<double> &mass_arr,  double &h);

void momentumEquation( vector<vector<unsigned>> &neighbours_matrix,  vector<double> &mass,  vector<vector<double>> &gradW_matrix, 
                      vector<vector<double>> &artificial_visc,  vector<double> &rho_arr, vector<double> p_arr,  string &state_equation_chosen);