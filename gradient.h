#include <stdio.h>
#include <vector>
#include "sorted_list.h"
using namespace std;

void gradW(const vector<double> &part_pos, const vector<vector<unsigned>> &neighbours_matrix, vector<vector<double>> &gradW_matrix, 
           double L[3], const double &h, const int &Nx, const int &Ny, const int &Nz);

double setArtificialViscosity(vector<vector<double>> &artificial_visc, const vector<double> &part_pos, const vector<vector<unsigned>> &neighbours_matrix, vector<double> &u_arr, 
                            const double &c_ab, const double &rho_ab, const double &alpha, const double &beta, const double &h);

void continuityEquation(const vector<double> &part_pos, const vector<double> &u_arr, const vector<vector<unsigned>> &neighbours_matrix, 
                        const vector<vector<double>> &gradW_matrix, vector<double> &drhodt_arr, vector<double> &rho_arr, const double &mass, const double &h);

void momentumEquation(const vector<vector<unsigned>> &neighbours_matrix, const vector<double> &mass, const vector<vector<double>> &gradW_matrix, 
                      vector<vector<double>> &artificial_visc, const vector<double> &rho_arr, vector<double> p_arr, const string &state_equation_chosen);