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


void findNeighbours(vector<double> &part_pos, vector<vector<unsigned>> &cell_pos,
                 vector<vector<unsigned>> &neighbours_matrix, vector<double> &L_d, const unsigned &Nx, const unsigned &Ny, const unsigned &Nz, const double &h, const int &kappa){

    // Sort all particles in their corresponding cell
    for (unsigned pos = 0; pos < part_pos.size()/3; pos ++){

        unsigned idx_i = part_pos[3*pos + 0] / (L_d[0] / Nx);
        unsigned idx_j = part_pos[3*pos + 1] / (L_d[1] / Ny);
        unsigned idx_k = part_pos[3*pos + 2] / (L_d[2] / Nz);
        
        idx_i = (idx_i == Nx) ? idx_i - 1 : idx_i;
        idx_j = (idx_j == Ny) ? idx_j - 1 : idx_j;
        idx_k = (idx_k == Nz) ? idx_k - 1 : idx_k;
        cell_pos[idx_i + Nx*idx_j + Ny*Nx*idx_k].push_back(pos);
    }
    

    // Find neighbours for each particle
    for (unsigned pos = 0; pos < part_pos.size()/3; pos ++){

        // Determine in which cell the particle is
        unsigned i_cell = part_pos[3*pos + 0] / (L_d[0] / Nx);
        unsigned j_cell = part_pos[3*pos + 1] / (L_d[1] / Ny);
        unsigned k_cell = part_pos[3*pos + 2] / (L_d[2] / Nz);

        i_cell = (i_cell >= Nx) ? Nx-1 : i_cell;
        j_cell = (j_cell >= Ny) ? Ny-1 : j_cell;
        k_cell = (k_cell >= Nx) ? Nz-1 : k_cell;

        // Define neighbouring cell indices
        unsigned i_inf = (i_cell == 0) ? 0 : i_cell - 1;
        unsigned i_supp = (i_cell < Nx - 1) ? i_cell +1 : (i_cell == Nx - 1) ? i_cell : i_cell - 1;

        unsigned j_inf = (j_cell == 0) ? 0 : j_cell - 1;
        unsigned j_supp = (j_cell < Ny - 1) ? j_cell +1 : (j_cell == Ny - 1) ? j_cell : j_cell - 1;

        unsigned k_inf = (k_cell == 0) ? 0 : k_cell - 1;
        unsigned k_supp = (k_cell < Nz - 1) ? k_cell +1 : (k_cell == Nz - 1) ? k_cell : k_cell - 1;

        // Iterate over (max) 26 adjacents cells to find neighbours 
        for (unsigned i = i_inf; i <= i_supp; i++){
            for (unsigned j = j_inf; j <= j_supp; j++){
                for (unsigned k = k_inf; k <= k_supp; k++){

                    vector<unsigned> &actual_cell = cell_pos[i + j*Nx + k*Nx*Ny];  

                    for (size_t idx_neighbour_it = 0 ; idx_neighbour_it < actual_cell.size(); idx_neighbour_it++) {
                        unsigned actual_cell_value = actual_cell[idx_neighbour_it]; 

                        if(actual_cell_value != pos){
                            
                            double rx, ry, rz, r2;
                            rx = (part_pos[3*pos] - part_pos[3*actual_cell_value])*(part_pos[3*pos] - part_pos[3*actual_cell_value]);
                            ry = (part_pos[3*pos + 1] - part_pos[3*actual_cell_value+1])*(part_pos[3*pos + 1] - part_pos[3*actual_cell_value+1]);
                            rz = (part_pos[3*pos + 2] - part_pos[3*actual_cell_value+2])*(part_pos[3*pos + 2] - part_pos[3*actual_cell_value+2]);
                            r2 = rx + ry + rz;
                            if(r2 <= kappa*kappa*h*h){
                                neighbours_matrix[pos].push_back(actual_cell_value); 
                             }
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
            ry = (part_pos[3*i+1] - part_pos[3*j+1] )*(part_pos[3*i+1] - part_pos[3*j+1] );
            rz = (part_pos[3*i+2] - part_pos[3*j+2] )*(part_pos[3*i+2] - part_pos[3*j+2] );
            r2 = rx + ry + rz;
            if(r2<= kappa*kappa*h*h){
                neighbours_matrix[i].push_back(j); 
                neighbours_matrix[j].push_back(i);
            }
        }
    }
}

void printNeighbours(vector<vector<unsigned>> &neighbours_matrix_1, vector<vector<unsigned>> &neighbours_matrix_2){

    for (unsigned i = 0; i < neighbours_matrix_1.size(); i++) {
        std::cout << "Particle " << i << " : ";

        std::cout << "{";
        for (unsigned j = 0; j < neighbours_matrix_1[i].size(); j++) {
            std::cout << neighbours_matrix_1[i][j];
            if (j != neighbours_matrix_1[i].size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "} (Linked-list) VS {";
        
        for (unsigned j = 0; j < neighbours_matrix_2[i].size(); j++) {
            std::cout << neighbours_matrix_2[i][j];
            if (j != neighbours_matrix_2[i].size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "} (naive)\n \n";
        
    }
}
