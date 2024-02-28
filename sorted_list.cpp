#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
#include <cstdlib>
#include <random>
#include <list>
#include <unordered_map>
#include <iostream>
using namespace std;


void linkedListAlgo(vector<double> &particle_x, vector<double> &particle_y, vector<double> &particle_z, vector<unsigned> &particle_i, vector<unsigned> &particle_j, vector<unsigned> &particle_k, 
                 vector<vector<int>> &neighbours_matrix, double L[3], const int &Nx, const int &Ny, const int &Nz, const double &h, const int &kappa){

    // Sort all particles in their corresponding cell
    for (unsigned ite = 0; ite < particle_x.size(); ite ++){
                int idx_i = particle_x[ite]/(L[0]/Nx);
                int idx_j = particle_y[ite]/(L[1]/Ny);
                int idx_k = particle_z[ite]/(L[2]/Nz);
                if (idx_i == Nx){
                    idx_i--;
                }
                if (idx_j == Ny){
                    idx_j--;
                }
                if (idx_k == Nz){
                    idx_k--;
                }
                particle_i.push_back(idx_i);
                particle_j.push_back(idx_j);
                particle_k.push_back(idx_k);
                
    }
    
    // Find neighbours for each particle
    for (unsigned x = 0; x < particle_x.size(); x++){
        
        // Determine in which cell the particle is a changer peut-être
        unsigned i_cell = particle_i[x];
        unsigned j_cell = particle_j[x];
        unsigned k_cell = particle_k[x];

        unsigned i_inf = i_cell -1;
        unsigned i_supp = i_cell+1;
        if(i_cell == 0 ){
                i_inf = 0;
                
        }
        
        if(i_cell == Nx-1 ){
                i_supp = i_cell;
                
                
                
        }

        unsigned j_inf = j_cell - 1;
        unsigned j_supp = j_cell + 1;
        if(j_cell == 0 ){
                j_inf = 0;
        }
        if(j_cell == Ny-1 ){
                j_supp = j_cell;
        }

        unsigned k_inf = k_cell -1;
        unsigned k_supp = k_cell+1;
        if(k_cell == 0 ){
                k_inf = 0;
        }
        if(k_cell == Nz-1 ){
                k_supp = k_cell;
        }

        // Iterate over all 26 adjacents cells to find neighbours 
        for (unsigned i = i_inf; i <= i_supp; i++){
            for (unsigned j = j_inf; j <= j_supp; j++){
                for (unsigned k = k_inf; k <= k_supp; k++){
                   
                    // Iterate over all particles to find the corresponding neighbours
                    for (unsigned l = x+1; l < particle_i.size(); l++){                         
                        if (particle_i[l] == i){
                            if (particle_j[l] == j){
                                if (particle_k[l] == k){
                                    double rx, ry, rz, r2;
                                    rx = (particle_x[x] - particle_x[l] )*(particle_x[x] - particle_x[l] );
                                    ry = (particle_y[x] - particle_y[l] )*(particle_y[x] - particle_y[l] );
                                    rz = (particle_z[x] - particle_z[l] )*(particle_z[x] - particle_z[l] );
                                    r2 = rx + ry + rz;
                                    
                                    if(r2<= kappa*kappa*h*h){
                                        neighbours_matrix[x].push_back(l); 
                                        neighbours_matrix[l].push_back(x);
                                    }            
                                }      
                            }
                        }        
                    }
                }
            }
       }
            
    }
} 

void naiveAlgo(vector<double> &particle_x, vector<double> &particle_y, vector<double> &particle_z, 
                 vector<vector<int>> &neighbours_matrix, const double &h, const int &kappa){

    // Find neighbours for each particle
    for (unsigned i = 0; i < particle_x.size(); i++){
        for (unsigned j = i+1; j < particle_x.size(); j++){

            double rx, ry, rz, r2;
            rx = (particle_x[i] - particle_x[j] )*(particle_x[i] - particle_x[j] );
            ry = (particle_y[i] - particle_y[j] )*(particle_y[i] - particle_y[j] );
            rz = (particle_z[i] - particle_z[j] )*(particle_z[i] - particle_z[j] );
            r2 = rx + ry + rz;
            if(r2<= kappa*kappa*h*h){
                neighbours_matrix[i].push_back(j); 
                neighbours_matrix[j].push_back(i);
            }
        }

    }
}
