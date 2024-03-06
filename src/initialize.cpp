#include <stdio.h>
#include <vector>
using namespace std;


void initializeMass(vector<double> &rho_arr,double s,vector<double> &mass_arr){

    double V = s*s*s;
    for(size_t i = 0; i < rho_arr.size(); i++){
        mass_arr[i] = rho_arr[i]*V;
    }    
}

void initializeRho(vector<double> &rho_arr,double &rho){

    for(size_t i = 0; i < rho_arr.size(); i++){
        rho_arr[i] = rho;
    }    
}

void initializeVelocity(vector<double> &u_arr,vector<double> &u_init){

    for(size_t i = 0; i < u_arr.size()/3; i++){
        u_arr[3*i] = u_init[0];
        u_arr[3*i+1] = u_init[1];
        u_arr[3*i+2] = u_init[2];
    }    
}
