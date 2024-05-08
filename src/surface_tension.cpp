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
                    vector<int> neighbours,
                    vector<vector<double>> gradW_matrix,
                    vector<double> mass,
                    vector<double> rho,
                    vector<double> pos,
                    vector<double> &F_vol){
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
            cout<<"dz = "<< dz << en    dl;

        }
    }
    */
    vector<double> normal(3*simParams.nb_moving_part,0.0);
    #pragma omp parallel for 
    for(int n = 0; n<simParams.nb_moving_part; n++){

        vector<double> &gradW = gradW_matrix[n];
        int size_neighbours = nb_neighbours[n];
        
        for(int idx = 0; idx<size_neighbours; idx++){ 

            int i_neig = neighbours[100*n + idx]; 
            double m_j = mass[i_neig];
            double rho_j = rho[i_neig];
            
            for( int coord = 0; coord <3; coord ++){

                double grad = gradW[3*idx+coord];
                normal[3*n+coord] += geomParams.h*m_j*grad/rho_j;
            }
        }
    }

    double alpha = simParams.alpha_st;
    #pragma omp parallel for 
    for(int n = 0; n<simParams.nb_moving_part; n++){

        int size_neighbours = nb_neighbours[n];
        
        for(int idx = 0; idx<size_neighbours; idx++){

            int i_neig = neighbours[100*n + idx];
            double K_ij = 2*thermoParam.rho_0/(rho[n]+rho[i_neig]);
            double r_ab = 0;
            vector<double> d_xyz(3);

            for (int coord = 0; coord < 3; coord++){
                
                d_xyz[coord] = pos[3 * n + coord] - pos[3 * i_neig + coord];
                r_ab += d_xyz[coord]*d_xyz[coord];
            }

            r_ab = sqrt(r_ab);
            double W_ab = W_coh(r_ab,geomParams.h);
            double m_a = mass[n];
            double m_b = mass[i_neig];
            double F_res = 0;

            for (int coord = 0; coord < 3; coord++){
                F_vol[3*n + coord] += -K_ij*(alpha * m_a * m_b * d_xyz[coord]*W_ab/r_ab 
                                  + alpha*(normal[3*n+coord]-normal[3*i_neig+coord]));
                F_res += F_vol[3*n + coord]*F_vol[3*n + coord];
            }

            simParams.F_st_max = sqrt(F_res);
        }
    }














    
}