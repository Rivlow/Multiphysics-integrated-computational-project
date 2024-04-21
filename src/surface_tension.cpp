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



void surfaceTension(SimulationData& simParams,
GeomData &geomParams,
vector<double> nb_neighbours,
vector<vector<int>> neighbours_matrix,
vector<double> mass,
vector<double> rho,
vector<double> pos,
vector<double> &F_vol
){
    
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
            /*cout<<"W = "<< W << endl;
            cout<<"dx = "<< dx << endl;
            cout<<"dy = "<< dy << endl;
            cout<<"dz = "<< dz << endl;*/

        }
    }
    
}