#include <stdio.h>
#include <vector>
#include <omp.h>

#include "gradient.h"
#include "find_neighbours.h"
#include "Kernel_functions.h"
#include "tools.h"
#include "structure.h"
#include "surface_tension.h"

using namespace std;

void surfaceTension(SimulationData& params,
vector<vector<double>> gradW_matrix,
vector<vector<int>> neighbours_matrix,
vector<double> mass,
vector<double> rho,
vector<double> pos,
vector<double> &F_vol
){
    int a = 2, b = 4;
    int nb_part = rho.size();
    vector<double> fprime(nb_part,1);
    vector<double> normal(3*nb_part,0);
    vector<double> div_normal(nb_part,0); 
    vector<double> div_star_normal(nb_part,0); 
    vector<double> N(nb_part,1);
    vector<double> C(nb_part,1);
    double h = params.h;
    double sigma = 72.8;

    #pragma omp parallel for
    for(int n=0; n<params.nb_moving_part ; n++){

        vector<double> &gradW = gradW_matrix[n];
        vector<int> &neighbours = neighbours_matrix[n];
        int size_neigbours = neighbours.size();
        double f = 0 ;
        for(int idx =0; idx<size_neigbours; idx++){
            double r2 = 0;
            for(int coord =0; coord < 3; coord++){
                double dcoord = pos[3*n+ coord] - pos[3*neighbours[idx]+coord];
                r2 = r2 + dcoord*dcoord;
            }
            double r = sqrt(r2);
            double W = f_cubic_spline(r,h);
            //cout<< "W = " << W << endl;
           // cout << "mass[neighbours[idx]]/rho[neighbours[idx]] = "<< mass[neighbours[idx]]/rho[neighbours[idx]] << endl;
            f = f + mass[neighbours[idx]]/rho[neighbours[idx]]*W;
            //cout<< " f = " << f << endl;
            
        }
        for(int exp = 0; exp < b; exp++){
                //cout<<"(a*f - a + 1) = "<< fprime[n] << endl;
                fprime[n] *= (a*f - a + 1);
            
            }
    }
    //printArray(fprime, fprime.size(),"fprime");
    #pragma omp parallel for
    for(int n=0; n<params.nb_moving_part ; n++){
        vector<double> &gradW = gradW_matrix[n];
        vector<int> &neighbours = neighbours_matrix[n];
        int size_neigbours = neighbours.size();
        //printArray(fprime, fprime.size(),"fprime");
        for(int idx = 0; idx < size_neigbours; idx++){
            double v = mass[neighbours[idx]]/rho[neighbours[idx]];
            double diff_fprime = fprime[neighbours[idx]] - fprime[n];
            //cout << " diif fprime " << diff_fprime << endl;
            //cout <<"particules_" << n;
            for(int coord = 0; coord < 3 ; coord++){
                 normal[3*n+coord] += v*diff_fprime*gradW[3*idx+coord];
                 //cout  << " coord_" << coord << " normal : " << normal[3*n+coord] << " et gradW " << gradW[3*idx+coord] << endl;
            }
            
        }
        
        double norm = 0;
        for(int coord = 0; coord<3; coord++){
            double d = normal[3*n+coord];
            norm += d*d;
            //cout << "norme : " << norm << endl;
        }
        if(norm >= 0.01/h){
            N[n] = 1;
        }
        else{
            N[n] = 0;
        }

    }
    //printArray(N, N.size(), "N");
    //printArray(normal, normal.size(),"normal");

    for(int n=0; n<params.nb_moving_part ; n++){
        vector<double> &gradW = gradW_matrix[n];
        vector<int> &neighbours = neighbours_matrix[n];
        int size_neigbours = neighbours.size();
        for(int idx = 0; idx<size_neigbours; idx++){
            double r2 = 0;
            for(int coord =0; coord < 3; coord++){
                double dcoord = pos[3*n+ coord] - pos[3*neighbours[idx]+coord];
                r2 = r2 + dcoord*dcoord;
            }
            double r = sqrt(r2);
            double W = f_cubic_spline(r,h);

            double v = mass[neighbours[idx]]/rho[neighbours[idx]];
            C[n] +=  v*N[n]*N[neighbours[idx]]*W;
        }
    }
    //printArray(C, C.size(), "C");

    #pragma omp parallel for
    for(int n=0; n<params.nb_moving_part ; n++){
        vector<double> &gradW = gradW_matrix[n];
        vector<int> &neighbours = neighbours_matrix[n];
        int size_neigbours = neighbours.size();
        for(int idx = 0; idx < size_neigbours; idx++){
            double dot = 0;
            for(int coord = 0; coord <3; coord++){
                dot += (normal[3*neighbours[idx]+coord] - normal[3*n+coord])*gradW[3*idx+coord];
                //cout << " grad_ "<< coord << " = "<< gradW[3*idx+coord] << endl; 
                //cout << " normal[3*neighbours[idx]+]_ "<< coord << " = "<< normal[3*neighbours[idx]+coord] << endl;
                //cout << "index " << idx << endl;
                //cout << " normal[3*n+]]_ "<< coord << " = "<< normal[3*n+coord] << endl; 
            }      
            //cout << " dot " << dot << endl;
            double v = mass[neighbours[idx]]/rho[neighbours[idx]];
            div_normal[n] += v*N[n]*N[neighbours[idx]]*dot;
        }
    div_normal[n] = div_normal[n]/C[n];
    
    for(int idx = 0; idx <size_neigbours; idx++){
        double r2 = 0;
        for(int coord =0; coord < 3; coord++){
            double dcoord = pos[3*n+ coord] - pos[3*neighbours[idx]+coord];
            r2 = r2 + dcoord*dcoord;
        }
        double r = sqrt(r2);
        double W = f_cubic_spline(r,h);
        double v = mass[neighbours[idx]]/rho[neighbours[idx]];
        div_star_normal[n] += v*W;

    }
    
    div_star_normal[n] = div_star_normal[n]*div_normal[n]/C[n];
    for(int coord = 0; coord <3; coord++){
        F_vol[3*n+coord]  = sigma/rho[n] *div_star_normal[n]*normal[3*n+coord]/C[n];
    } 
}








}