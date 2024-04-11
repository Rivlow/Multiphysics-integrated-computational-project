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

#include <omp.h>

using namespace std;

void sorted_list(SimulationData& params, 
                 vector<vector<int>> &cell_matrix,
                 vector<vector<int>> &neighbours_matrix,
                 vector<vector<double>> &gradW_matrix,
                 vector<vector<double>> &artificial_visc_matrix,
                 vector<double> &nb_neighbours,
                 vector<double> &pos){

    int Nx = params.Nx;
    int Ny = params.Ny;
    int Nz = params.Nz;
    int size_pos = pos.size() / 3;

    // Sort all particles in their corresponding cell
    for (int n = 0; n < size_pos; n++){

        int i = pos[3 * n + 0] / (params.L_d[0] / Nx);
        int j = pos[3 * n + 1] / (params.L_d[1] / Ny);
        int k = pos[3 * n + 2] / (params.L_d[2] / Nz);

        if (i < 0 || j < 0 || k < 0 || i > Nx || j > Ny || k > Ny){
            continue;
        }
        
        i = (i == Nx) ? i - 1 : i;
        j = (j == Ny) ? j - 1 : j;
        k = (k == Nz) ? k - 1 : k;
        cell_matrix[i + Nx * j + Ny * Nx * k].push_back(n);

    }

    // Find neighbours for each particle
    #pragma omp parallel for
    for (int n = 0; n < size_pos; n++){

        int i_cell = pos[3 * n + 0] / (params.L_d[0] / Nx);
        int j_cell = pos[3 * n + 1] / (params.L_d[1] / Ny);
        int k_cell = pos[3 * n + 2] / (params.L_d[2] / Nz);

        if (i_cell < 0 || j_cell < 0 || k_cell < 0 || i_cell > Nx || j_cell > Ny || k_cell > Ny){
            continue;
        }

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

        // Iterate over (max) 26 adjacents cells to find neighbours
        for (int i = i_inf; i <= i_supp; i++){
            for (int j = j_inf; j <= j_supp; j++){
                for (int k = k_inf; k <= k_supp; k++){

                    vector<int> &cell = cell_matrix[i + j * Nx + k * Nx * Ny];
                    int size_cell =cell.size();

                    // Iterate over neighbours
                    for (int idx = 0; idx < size_cell; idx++){

                        int idx_cell = cell[idx];

                        if (idx_cell != n){

                            double rx, ry, rz, r2;
                            rx = (pos[3 * n + 0] - pos[3 * idx_cell + 0]);
                            ry = (pos[3 * n + 1] - pos[3 * idx_cell + 1]);
                            rz = (pos[3 * n + 2] - pos[3 * idx_cell + 2]);
                            r2 = rx*rx + ry*ry + rz*rz;

                            int kappa = params.kappa;
                            double h = params.h;

                            if (r2 <= kappa * kappa *h * h){
    
                                neighbours_matrix[n].push_back(idx_cell);

                            }
                        }
                    }
                    gradW_matrix[n].resize(3*neighbours_matrix[n].size());
                    artificial_visc_matrix[n].resize(neighbours_matrix[n].size());
                    nb_neighbours[n] = neighbours_matrix[n].size();

                }
            }
        }
    }

    if (params.PRINT){
        cout << "findNeighbours passed" << endl;
    }

}

void naiveAlgo(SimulationData& params, 
               vector<vector<int>> &neighbours_matrix,
               vector<double> &pos){

    int nb_moving_part = params.nb_moving_part;

    // added by RB
    for (int i = 0; i < nb_moving_part; i++)
        neighbours_matrix[i].resize(0);


    // Find neighbours for each particle
    for (int i = 0; i < nb_moving_part; i++){

        for (int j = i + 1; j < nb_moving_part; j++){

            double rx, ry, rz, r2;
            rx = (pos[3 * i] - pos[3 * j]);
            ry = (pos[3 * i + 1] - pos[3 * j + 1]);
            rz = (pos[3 * i + 2] - pos[3 * j + 2]);
            r2 = rx*rx + ry*ry + rz*rz;

            int kappa = params.kappa;
            double h = params.h;

            if (r2 <= kappa * kappa *h * h){

                neighbours_matrix[i].push_back(j);
                neighbours_matrix[j].push_back(i);
            }
        }
    }
}

void printNeighbours(vector<vector<int>> &neighbours_matrix_linked,
                     vector<vector<int>> &neighbours_matrix_naive){

    for (int i = 0; i < int(neighbours_matrix_linked.size()); i++){
        std::cout << "Particle " << i << " : ";

        std::cout << "{";
        for (int j = 0; j < int(neighbours_matrix_linked[i].size()); j++){

            std::cout << neighbours_matrix_linked[i][j];
            if (j != int(neighbours_matrix_linked[i].size() - 1)){
                std::cout << ", ";
            }
        }
        std::cout << "} (Linked-list) VS {";

        for (int j = 0; j < int(neighbours_matrix_naive[i].size()); j++){

            std::cout << neighbours_matrix_naive[i][j];
            if (j != int(neighbours_matrix_naive[i].size() - 1)){
                std::cout << ", ";
            }
        }
        std::cout << "} (naive)\n \n";
    }
}

void CompareNeighbours( std::vector<std::vector<int>> &neighbours_matrix_linked,
                      std::vector<std::vector<int>> &neighbours_matrix_naive){
                        
    for (int i = 0; i < int(neighbours_matrix_linked.size()); i++){
        std::cout << "Particle " << i << " : ";

        std::cout << "{";
        for (int j = 0; j < int(neighbours_matrix_linked[i].size()); j++) {
            if (j < int(neighbours_matrix_naive[i].size())) {
                if (neighbours_matrix_linked[i][j] == neighbours_matrix_naive[i][j]){
                    std::cout << neighbours_matrix_linked[i][j];
                } else {
                    std::cout << "[" << neighbours_matrix_linked[i][j] << "-" << neighbours_matrix_naive[i][j] << "]";
                }
            } else {
                std::cout << "[" << neighbours_matrix_linked[i][j] << "-N/A]";
            }

            if (j != int(neighbours_matrix_linked[i].size() - 1)) {
                std::cout << ", ";
            }
        }
        std::cout << "} (Linked-list) VS {";

        for (int j = 0; j < int(neighbours_matrix_naive[i].size()); j++) {
            if (j >= int(neighbours_matrix_linked[i].size())) {
                std::cout << "[" << neighbours_matrix_naive[i][j] << "-N/A]";
            }

            if (j != int(neighbours_matrix_naive[i].size() - 1)) {
                std::cout << ", ";
            }
        }
        std::cout << "} (naive)\n \n";
    }
}
