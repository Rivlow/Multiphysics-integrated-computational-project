#include <stdio.h>
#include <vector>
#include <cmath>
#include <omp.h>

#include "gradient.h"
#include "find_neighbours.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"
#include "surface_tension.h"

using namespace std;

double W_coh(double r, double h){
    double W = 0.0;
    double cst = 32/(M_PI*h*h*h*h*h*h*h*h*h);
    if(2.0*r>h && r<=h){
        W = cst*(h-r)*(h-r)*(h-r)*r*r*r;
    }
    else if(r>0 && 2.0*r<=h){
        W = cst *( 2.0*(h-r)*(h-r)*(h-r)*r*r*r - h*h*h*h*h*h/64);
    }
    else{
        W = 0.0;
    }
    return W;
}

void surfaceTension(SimulationData& simParams,
GeomData &geomParams,
ThermoData &thermoParam,
vector<double> nb_neighbours,
vector<vector<int>> neighbours_matrix,
vector<vector<double>> gradW_matrix,
vector<double> mass,
vector<double> rho,
vector<double> pos,
vector<double> &F_vol
){
    /*
    double alpha = 12.5;
    #pragma omp parallel for 
    for(int n = 0; n<simParams.nb_moving_part; n++){
        vector<int> &neighbours = neighbours_matrix[n];
        int size_neighbours = nb_neighbours[n];
        
        for(int idx = 0; idx<size_neighbours; idx++)
        {
            
            double dx = pos[3*n+0] - pos[3*neighbours[idx]+0];
            double dy = pos[3*n+1] - pos[3*neighbours[idx]+1];
            double dz = pos[3*n+2] - pos[3*neighbours[idx]+2];
            double r2 = dx*dx + dy*dy + dz*dz;
            double r = sqrt(r2);
            double W = f_cubic_spline(r,geomParams.h);
            
            F_vol[3*n + 0] += -(alpha/mass[n])*mass[neighbours[idx]]*dx*W;
            F_vol[3*n + 1] += -(alpha/mass[n])*mass[neighbours[idx]]*dy*W;
            F_vol[3*n + 2] += -(alpha/mass[n])*mass[neighbours[idx]]*dz*W;
            cout<<"W = "<< W << endl;
            cout<<"dx = "<< dx << endl;
            cout<<"dy = "<< dy << endl;
            cout<<"dz = "<< dz << endl;

        }
    }
    */
    vector<double> normal(3*simParams.nb_moving_part,0.0);
    #pragma omp parallel for 
    for(int n = 0; n<simParams.nb_moving_part; n++){
        vector<int> &neighbours = neighbours_matrix[n];
        vector<double> &gradW = gradW_matrix[n];
        int size_neighbours = nb_neighbours[n];
        
        for(int idx = 0; idx<size_neighbours; idx++)
        {         
            int i_neig = neighbours[idx]; 
            double m_j = mass[i_neig];
            double rho_j = rho[i_neig];
            
            for( int coord = 0; coord <3; coord ++){
                double grad = gradW[3*idx+coord];
                normal[3*n+coord] += geomParams.h*m_j*grad/rho_j;

            }
        }
    }


    double alpha = 100;
    #pragma omp parallel for 
    for(int n = 0; n<simParams.nb_moving_part; n++){
        vector<int> &neighbours = neighbours_matrix[n];
        int size_neighbours = nb_neighbours[n];
        
        for(int idx = 0; idx<size_neighbours; idx++)
        {
            double K_ij = 2*thermoParam.rho_0/(rho[n]+rho[neighbours[idx]]);
            double dx = pos[3*n+0] - pos[3*neighbours[idx]+0];
            double dy = pos[3*n+1] - pos[3*neighbours[idx]+1];
            double dz = pos[3*n+2] - pos[3*neighbours[idx]+2];
            double r2 = dx*dx + dy*dy + dz*dz;
            double r = sqrt(r2);
            double W = W_coh(r,geomParams.h);

            
            F_vol[3*n + 0] += -K_ij*((alpha*mass[n])*mass[neighbours[idx]]*dx*W/r + alpha*(normal[3*n+0]-normal[3*neighbours[idx]+0]));
            F_vol[3*n + 1] += -K_ij*((alpha*mass[n])*mass[neighbours[idx]]*dy*W/r + alpha*(normal[3*n+1]-normal[3*neighbours[idx]+1]));
            F_vol[3*n + 2] += -K_ij*((alpha*mass[n])*mass[neighbours[idx]]*dz*W/r + alpha*(normal[3*n+2]-normal[3*neighbours[idx]+2]));
            
             
        }
    }














    
}