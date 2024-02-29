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


void sort_tableaux(vector<unsigned>& tableau_1, vector<double>& tableau_2) {
    // Créer un vecteur de paires pour associer chaque élément de tableau_1 avec son indice original
    vector<pair<int, int>> temp;
    for (int i = 0; i < tableau_1.size(); ++i) {
        temp.push_back(make_pair(tableau_1[i], i));
    }

    // Trier temp en fonction de la première valeur de chaque paire (les valeurs de tableau_1)
    sort(temp.begin(), temp.end());

    // Réorganiser tableau_1 en fonction de l'ordre obtenu et réorganiser tableau_2 en même temps
    vector<double> temp_tableau_2(tableau_2.size());
    for (int i = 0; i < tableau_1.size(); ++i) {
        tableau_1[i] = temp[i].first;
        int index_tableau_1 = temp[i].second;
        temp_tableau_2[i * 3] = tableau_2[index_tableau_1 * 3];
        temp_tableau_2[i * 3 + 1] = tableau_2[index_tableau_1 * 3 + 1];
        temp_tableau_2[i * 3 + 2] = tableau_2[index_tableau_1 * 3 + 2];
    }

    // Copier les valeurs réorganisées de temp_tableau_2 dans tableau_2
    tableau_2 = temp_tableau_2;
}



void linkedListAlgo(vector<double> &part_pos, vector<unsigned> &cell_pos, vector<unsigned> &tab_cumul,
                 vector<vector<int>> &neighbours_matrix, double L[3], const int &Nx, const int &Ny, const int &Nz, const double &h, const int &kappa){

    // Sort all particles in their corresponding cell
    for (unsigned ite = 0; ite < part_pos.size()/3; ite ++){

        int idx_i = part_pos[3*ite + 0] / (L[0] / Nx);
        int idx_j = part_pos[3*ite + 1] / (L[1] / Ny);
        int idx_k = part_pos[3*ite + 2] / (L[2] / Nz);

        idx_i = (idx_i == Nx) ? idx_i - 1 : idx_i;
        idx_j = (idx_j == Ny) ? idx_j - 1 : idx_j;
        idx_k = (idx_k == Nz) ? idx_k - 1 : idx_k;

        cell_pos.push_back(idx_i + Nx*idx_j + Ny*Nx*idx_k);
        printf("%d \n", idx_i + Nx*idx_j + Ny*Nx*idx_k);

        for (unsigned k = Nx*Ny*Nz - 1; k > idx_i + Nx*idx_j + Ny*Nx*idx_k; k--){
            tab_cumul[k]++;
        }
    }

    sort_tableaux(cell_pos, part_pos);

    // Find neighbours for each particle
    for (unsigned x = 0; x < part_pos.size()/3; x++){
        
        // Determine in which cell the particle is
        int a = cell_pos[x]%(Ny*Nx);
        unsigned i_cell = a%Nx;
        unsigned j_cell = (a - i_cell)/Nx;
        unsigned k_cell = (cell_pos[x] - a)/(Ny*Nx);

        // Define neighboring cell indices
        unsigned i_inf = (i_cell == 0) ? 0 : i_cell - 1;
        unsigned i_supp = (i_cell == Nx - 1) ? i_cell : i_cell + 1;

        unsigned j_inf = (j_cell == 0) ? 0 : j_cell - 1;
        unsigned j_supp = (j_cell == Ny - 1) ? j_cell : j_cell + 1;

        unsigned k_inf = (k_cell == 0) ? 0 : k_cell - 1;
        unsigned k_supp = (k_cell == Nz - 1) ? k_cell : k_cell + 1;
        
        // Iterate over all 26 adjacents cells to find neighbours 
        for (unsigned i = i_inf; i <= i_supp; i++){
            for (unsigned j = j_inf; j <= j_supp; j++){
                for (unsigned k = k_inf; k <= k_supp; k++){
                   
                    // Iterate over all particles to find the corresponding neighbours   
                    unsigned actual_cell = cell_pos[i + j*Nx + k*Nx*Ny];   
                    unsigned range_to_look = tab_cumul[actual_cell+1] - tab_cumul[actual_cell];
                    unsigned idx_init_to_look = cell_pos[tab_cumul[actual_cell]];
                    unsigned idx_end_to_look = cell_pos[tab_cumul[actual_cell+1]];

                    for (unsigned init = idx_init_to_look; init < idx_end_to_look, init++){
                        double rx, ry, rz, r2;
                        rx = (part_pos[3*x] - part_pos[l] )*(part_pos[x] - part_pos[l] );
                        ry = (part_pos[3*x + 1] - part_pos[l] )*(part_pos[x] - part_pos[l] );
                        rz = (part_pos[3*x + 2] - part_pos[l] )*(part_pos[x] - part_pos[l] );
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
