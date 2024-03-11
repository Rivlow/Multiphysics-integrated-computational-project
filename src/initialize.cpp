#include <stdio.h>
#include <vector>
#include <string.h>
#include <string>
using namespace std;


void initializeMass(vector<double> &rho_arr,double &s, vector<double> &mass_arr){

    double V = s*s*s;
    for(size_t i = 0; i < rho_arr.size(); i++){
        mass_arr[i] = rho_arr[i]*V;
    }    
}

void initializeRho(vector<double> &rho_arr,double &rho, double &rho_0, double &c_0, 
                   double &M, double &g, double &H, double &R, double &T, double &gamma, string &state_equation_chosen, 
                   string &state_initial_condition){

    if (state_initial_condition == "Hydrostatic"){

        if (state_equation_chosen == "Ideal gaz law"){

            for(size_t i = 0; i < rho_arr.size(); i++){
                rho_arr[i] = rho_0*(1+ M*rho_0*g*H/(R*T));
            }    
        }
        if (state_equation_chosen == "Quasi incompresible fluid"){

            double B = c_0*c_0*rho_0/gamma;
            for(size_t i = 0; i < rho_arr.size(); i++){
                
                rho_arr[i] = rho_0*(1+ rho_0*g*H/B);
            }    
        }
    }

    else{

        for(size_t i = 0; i < rho_arr.size(); i++){
        rho_arr[i] = rho;
        }    
    }

    
}

void initializeVelocity(vector<double> &u_arr,vector<double> &u_init){

    for(size_t i = 0; i < u_arr.size()/3; i++){
        u_arr[3*i] = u_init[0];
        u_arr[3*i+1] = u_init[1];
        u_arr[3*i+2] = u_init[2];
    }    
}
