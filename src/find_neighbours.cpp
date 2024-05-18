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

#include "find_neighbours.h"
#include "structure.h"
#include "tools.h"

#include <omp.h>

using namespace std;

void sortedList(GeomData &geomParams,
                SimulationData &simParams, 
                vector<vector<int>> &cell_matrix,
                vector<int> &neighbours,
                vector<int> &track_surface,
                vector<vector<double>> &gradW_matrix,
                vector<vector<double>> &W_matrix,
                vector<vector<double>> &pi_matrix,
                vector<double> &nb_neighbours,
                vector<double> &type,
                vector<double> &pos){
    
    int Nx = geomParams.Nx;
    int Ny = geomParams.Ny;
    int Nz = geomParams.Nz;
    double Lx = geomParams.L_d[0], Ly = geomParams.L_d[1], Lz = geomParams.L_d[2];
    int size_pos = pos.size() / 3;
    // Sort all particles in their corresponding cell
    for (int n = 0; n < size_pos; n++){

        int i = pos[3 * n + 0] / (Lx / Nx);
        int j = pos[3 * n + 1] / (Ly / Ny);
        int k = pos[3 * n + 2] / (Lz / Nz);
        
        // Skip particules outside of the domain
        if (i < 0 || j < 0 || k < 0 || i > Nx || j > Ny || k > Nz) continue;
        
        // Modify index of particules at boundaries
        i = (i == Nx) ? i - 1 : i;
        j = (j == Ny) ? j - 1 : j;
        k = (k == Nz) ? k - 1 : k;
        
        cell_matrix[i + Nx * j + Ny * Nx * k].push_back(n);
        
        
    }

    
    // Find neighbours for each particle
    #pragma omp parallel for
    for (int n = 0; n < size_pos; n++){

        int i_cell = pos[3 * n + 0] / (geomParams.L_d[0] / Nx);
        int j_cell = pos[3 * n + 1] / (geomParams.L_d[1] / Ny);
        int k_cell = pos[3 * n + 2] / (geomParams.L_d[2] / Nz);

        // Skip particules outside of the domain
        if (i_cell < 0 || j_cell < 0 || k_cell < 0 || i_cell > Nx || j_cell > Ny || k_cell > Nz) continue;
        

        // Modify index of particules at boundaries
        i_cell = (i_cell >= Nx) ? Nx - 1 : i_cell;
        j_cell = (j_cell >= Ny) ? Ny - 1 : j_cell;
        k_cell = (k_cell >= Nx) ? Nz - 1 : k_cell;

        // Define neighbouring cell indices
        int i_inf = (i_cell == 0) ? 0 : i_cell - 1;
        int i_supp = (i_cell < Nx - 1) ? i_cell + 1 : (i_cell == Nx - 1) ? i_cell
                                                                         : i_cell - 1;
        int j_inf = (j_cell == 0) ? 0 : j_cell - 1;
        int j_supp = (j_cell < Ny - 1) ? j_cell + 1 : (j_cell == Ny - 1) ? j_cell
                                                                         : j_cell - 1;
        int k_inf = (k_cell == 0) ? 0 : k_cell - 1;
        int k_supp = (k_cell < Nz - 1) ? k_cell + 1 : (k_cell == Nz - 1) ? k_cell
                                                                         : k_cell - 1;
        int it = 0;

        // Iterate over (max) 26 adjacents cells to find neighbours
        for (int k = k_inf; k <= k_supp; k++){
            for (int j = j_inf; j <= j_supp; j++){
                for (int i = i_inf; i <= i_supp; i++){

                    vector<int> &cell = cell_matrix[i + j * Nx + k * Nx * Ny];
                    int size_cell = cell.size();

                    // Iterate over particles in cell
                    for (int idx = 0; idx < size_cell; idx++){

                        int idx_cell = cell[idx];

                        if (type[idx_cell] == 2.0) continue;
                        else{

                            if (idx_cell != n){

                                double rx, ry, rz, r2;
                                rx = (pos[3 * n + 0] - pos[3 * idx_cell + 0]);
                                ry = (pos[3 * n + 1] - pos[3 * idx_cell + 1]);
                                rz = (pos[3 * n + 2] - pos[3 * idx_cell + 2]);
                                r2 = rx*rx + ry*ry + rz*rz;

                                int kappa = geomParams.kappa;
                                double h = geomParams.h;

                                if (r2 <= kappa * kappa *h * h){
                                    
                                    neighbours[100*n + it++] = idx_cell;

                                   /* if (type[n] == 1){ //FP
           

                                        double theta = atan2(rz, rx); 
                                        //double phi = atan2(ry, rx); 
                                        

                                        theta = theta * 180.0 / M_PI  + (theta < 0 ? 360.0 : 0);
                                        //phi = phi * 180.0 / M_PI + (phi < 0 ? 360.0 : 0) ;// phi between 0 et 360 degrees
                                       
                                        int i_idx = static_cast<int>(theta / 45.0);
                                        //int j_idx = static_cast<int>(phi / 45.0);
                                        
                                        
                                        track_surface[8*n + i_idx]++;
                                    }*/
                                    
                                    
                                }
                            
                            }
                        }
                    }
                }
            }
        }
        
        gradW_matrix[n].resize(3*it);
        W_matrix[n].resize(it);
        pi_matrix[n].resize(it);
        nb_neighbours[n] = it;
       
    }

    if (simParams.PRINT) cout << "findNeighbours passed" << endl;


}

/*
void naiveAlgo(GeomData &geomParams,
               SimulationData &simParams, 
               vector<vector<int>> &neighbours_matrix,
               vector<double> &pos){

    int pos_size = pos.size()/3; 

    // added by RB
    //for (int i = 0; i < nb_moving_part; i++)
        //neighbours_matrix[i].resize(0);

    int it_i = 0, it_j = 0;

    // Find neighbours for each particle
    for (int i = 0; i < pos_size; i++){
        for (int j = i + 1; j < pos_size; j++){

            double rx, ry, rz, r2;
            rx = (pos[3 * i] - pos[3 * j]);
            ry = (pos[3 * i + 1] - pos[3 * j + 1]);
            rz = (pos[3 * i + 2] - pos[3 * j + 2]);
            r2 = rx*rx + ry*ry + rz*rz;

            int kappa = geomParams.kappa;
            double h = geomParams.h;

            if (r2 <= kappa * kappa *h * h){
                neighbours_matrix[i][it_i++];
                neighbours_matrix[j][it_j++];
            }
        }
    }
}
*/

void printNeighbours(vector<vector<int>> &neighbours_matrix_linked,
                     vector<vector<int>> &neighbours_matrix_naive){

    for (int i = 0; i < int(neighbours_matrix_linked.size()); i++){
        cout << "Particle " << i << " : ";

        cout << "{";
        for (int j = 0; j < int(neighbours_matrix_linked[i].size()); j++){

            cout << neighbours_matrix_linked[i][j];
            if (j != int(neighbours_matrix_linked[i].size() - 1))
                cout << ", ";
            
        }
        cout << "} (Linked-list) VS {";

        for (int j = 0; j < int(neighbours_matrix_naive[i].size()); j++){

            cout << neighbours_matrix_naive[i][j];
            if (j != int(neighbours_matrix_naive[i].size() - 1))
                cout << ", ";
            
        }
        cout << "} (naive)\n \n";
    }
}

void CompareNeighbours( vector<vector<int>> &neighbours_matrix_linked,
                      vector<vector<int>> &neighbours_matrix_naive){
                        
    for (int i = 0; i < int(neighbours_matrix_linked.size()); i++){
        cout << "Particle " << i << " : ";

        cout << "{";
        for (int j = 0; j < int(neighbours_matrix_linked[i].size()); j++) {
            if (j < int(neighbours_matrix_naive[i].size())) {
                if (neighbours_matrix_linked[i][j] == neighbours_matrix_naive[i][j]){
                    cout << neighbours_matrix_linked[i][j];
                } else {
                    cout << "[" << neighbours_matrix_linked[i][j] << "-" << neighbours_matrix_naive[i][j] << "]";
                }
            } else {
                cout << "[" << neighbours_matrix_linked[i][j] << "-N/A]";
            }

            if (j != int(neighbours_matrix_linked[i].size() - 1)) {
                cout << ", ";
            }
        }
        cout << "} (Linked-list) VS {";

        for (int j = 0; j < int(neighbours_matrix_naive[i].size()); j++) {
            if (j >= int(neighbours_matrix_linked[i].size())) {
                cout << "[" << neighbours_matrix_naive[i][j] << "-N/A]";
            }

            if (j != int(neighbours_matrix_naive[i].size() - 1)) {
                cout << ", ";
            }
        }
        cout << "} (naive)\n \n";
    }
}
