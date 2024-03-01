#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
#include <cstdlib>
#include <random>
#include <list>
#include <unordered_map>
#include <iostream>
#include <chrono>
#include <algorithm>

using namespace std;


void linkedListAlgo(vector<double> &part_pos, vector<vector<unsigned>> &cell_pos,
                 vector<vector<unsigned>> &neighbours_matrix, double L[3], const unsigned &Nx, const unsigned &Ny, const unsigned &Nz, const double &h, const int &kappa){

    // Sort all particles in their corresponding cell
    for (unsigned pos = 0; pos < part_pos.size()/3; pos ++){

        unsigned idx_i = part_pos[3*pos + 0] / (L[0] / Nx);
        unsigned idx_j = part_pos[3*pos + 1] / (L[1] / Ny);
        unsigned idx_k = part_pos[3*pos + 2] / (L[2] / Nz);

        idx_i = (idx_i == Nx) ? idx_i - 1 : idx_i;
        idx_j = (idx_j == Ny) ? idx_j - 1 : idx_j;
        idx_k = (idx_k == Nz) ? idx_k - 1 : idx_k;

        cell_pos[idx_i + Nx*idx_j + Ny*Nx*idx_k].push_back(pos);
    }

    // Find neighbours for each particle
    for (unsigned pos = 0; pos < part_pos.size()/3; pos ++){
        
        // Determine in which cell the particle is
        unsigned i_cell = part_pos[3*pos + 0] / (L[0] / Nx);
        unsigned j_cell = part_pos[3*pos + 1] / (L[1] / Ny);
        unsigned k_cell = part_pos[3*pos + 2] / (L[2] / Nz);

        // Define neighbouring cell indices
        unsigned i_inf = (i_cell == 0) ? 0 : i_cell - 1;
        unsigned i_supp = (i_cell < Nx - 1) ? i_cell +1 : (i_cell == Nx - 1) ? i_cell : i_cell - 1;

        unsigned j_inf = (j_cell == 0) ? 0 : j_cell - 1;
        unsigned j_supp = (j_cell < Ny - 1) ? j_cell +1 : (j_cell == Ny - 1) ? j_cell : j_cell - 1;

        unsigned k_inf = (k_cell == 0) ? 0 : k_cell - 1;
        unsigned k_supp = (k_cell < Nz - 1) ? k_cell +1 : (k_cell == Nz - 1) ? k_cell : k_cell - 1;

        
        printf("actual cell : (%d, %d, %d) -> ", i_cell, j_cell, k_cell);
        printf("neighbour lower : (%d, %d, %d) ", i_inf, j_inf, k_inf);
        printf("neighbour upper : (%d, %d, %d) ", i_supp, j_supp, k_supp);
        
        // Iterate over (max) 26 adjacents cells to find neighbours 
        for (unsigned i = i_inf; i <= i_supp; i++){
            for (unsigned j = j_inf; j <= j_supp; j++){
                for (unsigned k = k_inf; k <= k_supp; k++){

                    vector<unsigned> &actual_cell = cell_pos[i + j*Nx + k*Nx*Ny];  
                    int val = i + j*Nx + k*Nx*Ny;
                    printf("cell en question : %d \n \n", val);

                    for (unsigned idx_neighbour = 0; idx_neighbour < actual_cell.size(); idx_neighbour++){
                        //printf("idx_neighbour = %d \n", idx_neighbour);

                        double rx, ry, rz, r2;
                        rx = (part_pos[3*pos] - part_pos[3*idx_neighbour])*(part_pos[3*pos] - part_pos[3*idx_neighbour]);
                        ry = (part_pos[3*pos + 1] - part_pos[3*idx_neighbour+1])*(part_pos[3*pos + 1] - part_pos[3*idx_neighbour+1]);
                        rz = (part_pos[3*pos + 2] - part_pos[3*idx_neighbour+2])*(part_pos[3*pos + 2] - part_pos[3*idx_neighbour+2]);
                        r2 = rx + ry + rz;
                        
                        if(r2<= kappa*kappa*h*h){
                            neighbours_matrix[pos].push_back(idx_neighbour); 
                            neighbours_matrix[idx_neighbour].push_back(pos);
                        }      
                    }              
                }
            }
        } 
    }
} 

void naiveAlgo(vector<double> &part_pos, vector<vector<unsigned>> &neighbours_matrix, const double &h, const int &kappa){

    // Find neighbours for each particle
    for (unsigned i = 0; i < part_pos.size()/3; i++){
        for (unsigned j = i+1; j < part_pos.size()/3; j++){


            double rx, ry, rz, r2;
            rx = (part_pos[3*i] - part_pos[3*j] )*(part_pos[3*i] - part_pos[3*j] );
            ry = (part_pos[3*i+1] - part_pos[3*j+1] )*(part_pos[3*i+1] - part_pos[3*j]+1 );
            rz = (part_pos[3*i+2] - part_pos[3*j+2] )*(part_pos[3*i+2] - part_pos[3*j+2] );
            r2 = rx + ry + rz;
            if(r2<= kappa*kappa*h*h){
                neighbours_matrix[i].push_back(j); 
                neighbours_matrix[j].push_back(i);
            }
        }
    }
}
