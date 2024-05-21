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
#include "Kernel.h"

using namespace std;


void surfaceTension(SimulationData& simParams,
                    GeomData &geomParams,
                    ThermoData &thermoParam,
                    vector<double> nb_neighbours,
                    vector<int> neighbours,
                    vector<int> &track_surface,
                    vector<double> &N_smoothed,
                    vector<vector<double>> gradW_matrix,
                    vector<vector<double>> W_matrix,
                    vector<double> mass,
                    vector<double> rho,
                    vector<double> pos,
                    vector<double> &F_vol,
                    vector<double> type){

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
            
            if(type[neighbours[100*n + idx]] == 1){
                int i_neig = neighbours[100*n + idx];
                double K_ij = 2*thermoParam.rho_0/(rho[n]+rho[i_neig]);
                double r_ab = 0;
                vector<double> d_xyz(3);
                
                for (int coord = 0; coord < 3; coord++){
                    
                    d_xyz[coord] = pos[3 * n + coord] - pos[3 * i_neig + coord];
                    r_ab += d_xyz[coord]*d_xyz[coord];
                }
                
                r_ab = sqrt(r_ab);
                double W_ab = W_coh(r_ab,geomParams.kappa*geomParams.h, simParams);
                
                double m_a = mass[n];
                double m_b = mass[i_neig];
                double F_res = 0;
                
                for (int coord = 0; coord < 3; coord++){
                    
                    F_vol[3*n + coord] += -K_ij*(alpha * m_a * m_b * d_xyz[coord]*W_ab/r_ab 
                                    + alpha*(normal[3*n+coord]-normal[3*i_neig+coord]));

                    //cout << "F_vol[3*n + coord] (in surface tension)" << -K_ij*(alpha * m_a * m_b * d_xyz[coord]*W_ab/r_ab 
                             //       + alpha*(normal[3*n+coord]-normal[3*i_neig+coord])) << endl;
                
                    F_res += F_vol[3*n + coord]*F_vol[3*n + coord];
                }
                
                simParams.F_st_max = sqrt(F_res);
            }
        }     
    }  
}
