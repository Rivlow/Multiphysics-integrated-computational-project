#include <stdio.h>
#include <vector>
#include <omp.h>

#include "gradient.h"
#include "find_neighbours.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"
#include "surface_tension.h"

using namespace std;

void surfaceParticle(SimulationData& params, vector<vector<int>> neighbours_matrix,vector<double> pos, vector<double> mass,
                     vector<double> rho, vector<vector<double>> gradW_matrix,vector<int> &surface_part){
    
    for(int n = 0; n < params.nb_moving_part ; n++){
        vector<int> &neighbours = neighbours_matrix[n];
        
        for(int idx = 0; idx < neighbours.size() ; idx++){

            double distz = pos[3*neighbours[idx] + 2] - pos[3*n + 2];
            if(distz >= 0.0){
            
                double distx = pos[3*neighbours[idx] + 0] - pos[3*n + 0];
                double disty = pos[3*neighbours[idx] + 1] - pos[3*n + 1];
            }
        }
    }
}

void surfaceTension(SimulationData& params,
vector<vector<double>> gradW_matrix,
vector<vector<int>> neighbours_matrix,
vector<double> mass,
vector<double> rho,
vector<double> pos,
vector<double> &F_vol
){

    
}