#include <stdio.h>
#include <vector>
#include <string>
#include "initialize.h"


using namespace std;


void initializeMass(vector<double> &rho_arr,double &s, vector<double> &mass_arr){

    double V = s*s*s;
    for(size_t i = 0; i < rho_arr.size(); i++){
        mass_arr[i] = rho_arr[i]*V;
    }    
}

void initializeRho(size_t nb_moving_part, vector<double> &pos_arr, vector<double> &rho_arr,double &rho_moving, double &rho_fixed, double &rho_0, double &c_0, 
                   double &M, double &g, double &R, double &T, double &gamma, std::string &state_equation_chosen, 
                   std::string &state_initial_condition){

    if (state_initial_condition == "Hydrostatic"){

        if (state_equation_chosen == "Ideal gaz law"){

            for(size_t i = 0; i < rho_arr.size(); i++){
                rho_arr[i] = (i < nb_moving_part) ? rho_0*(1+ M*rho_0*g*pos_arr[3*i + 2]/(R*T)) : rho_fixed;
            }    
        }
        if (state_equation_chosen == "Quasi incompresible fluid"){

            double B = c_0*c_0*rho_0/gamma;
            for(size_t i = 0; i < rho_arr.size(); i++){
                
                rho_arr[i] = (i < nb_moving_part) ? rho_0*(1+ rho_0*g*pos_arr[3*i + 2]/B) : rho_fixed;
            }    
        }
    }
    else{

        for(size_t i = 0; i < rho_arr.size(); i++){
        rho_arr[i] = (i < nb_moving_part) ? rho_moving : rho_fixed;
        }    
    }
}

void initializeVelocity(size_t nb_moving_part, vector<double> &u_arr,vector<double> &u_init){

    for(size_t i = 0; i < u_arr.size()/3; i++){

        if (i < nb_moving_part){
            u_arr[3*i] = u_init[0];
            u_arr[3*i+1] = u_init[1];
            u_arr[3*i+2] = u_init[2];

        }
        else{
            u_arr[3*i] = 0;
            u_arr[3*i+1] = 0;
            u_arr[3*i+2] = 0;
        }
    }    
}

void initializeViscosity(vector<vector<double>> &artificial_visc_matrix){

    for(size_t i = 0; i<artificial_visc_matrix.size();i++){
        for(size_t j=0 ; artificial_visc_matrix[i].size();j++){
            artificial_visc_matrix[i][j] = 0.0;
        }
    }  
}
