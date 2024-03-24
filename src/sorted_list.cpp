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
#include "sorted_list.h"

using namespace std;

void findNeighbours(size_t nb_moving_part, vector<double> &pos_arr,
 vector<vector<int>> &cell_matrix,
                    vector<vector<int>> &neighbours_matrix, 
                    vector<double> &L_d, const int &Nx, 
                    const int &Ny, const int &Nz,
                     const double &h, const int &kappa)
{

    // Sort all particles in their corresponding cell
    for (size_t pos = 0; pos < pos_arr.size() / 3; pos++)
    {

        int idx_i = pos_arr[3 * pos + 0] / (L_d[0] / Nx);
        int idx_j = pos_arr[3 * pos + 1] / (L_d[1] / Ny);
        int idx_k = pos_arr[3 * pos + 2] / (L_d[2] / Nz);

        idx_i = (idx_i == Nx) ? idx_i - 1 : idx_i;
        idx_j = (idx_j == Ny) ? idx_j - 1 : idx_j;
        idx_k = (idx_k == Nz) ? idx_k - 1 : idx_k;
        cell_matrix[idx_i + Nx * idx_j + Ny * Nx * idx_k].push_back(pos);

        // cout << "For part : " << pos << ", cell's index = (" << idx_i << ", " << idx_j << ", " << idx_k << ")" << endl;
    }

    // Find neighbours for each particle
    for (int pos = 0; pos < nb_moving_part; pos++)
    {

        // Determine in which cell the particle is
        int i_cell = pos_arr[3 * pos + 0] / (L_d[0] / Nx);
        int j_cell = pos_arr[3 * pos + 1] / (L_d[1] / Ny);
        int k_cell = pos_arr[3 * pos + 2] / (L_d[2] / Nz);

        // cout << "cell's indices computed" << endl;

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
        // cout << "cell's neighbours computed" << endl;

        // Iterate over (max) 26 adjacents cells to find neighbours
        for (size_t i = i_inf; i <= i_supp; i++)
        {
            for (size_t j = j_inf; j <= j_supp; j++)
            {
                for (size_t k = k_inf; k <= k_supp; k++)
                {
                    vector<int> &actual_cell = cell_matrix[i + j * Nx + k * Nx * Ny];
                    // cout << "len(actual_cell vector) : " << actual_cell.size() << endl;

                    if (actual_cell.size() > 0)
                    {
                        for (size_t idx_neighbour_it = 0; idx_neighbour_it < actual_cell.size(); idx_neighbour_it++)
                        {
                            int actual_cell_value = actual_cell[idx_neighbour_it];
                            // cout << "actual_cell_value defined" << endl;

                            if (actual_cell_value != pos)
                            {
                                double rx, ry, rz, r2;
                                rx = (pos_arr[3 * pos] - pos_arr[3 * actual_cell_value]) * (pos_arr[3 * pos] - pos_arr[3 * actual_cell_value]);
                                ry = (pos_arr[3 * pos + 1] - pos_arr[3 * actual_cell_value + 1]) * (pos_arr[3 * pos + 1] - pos_arr[3 * actual_cell_value + 1]);
                                rz = (pos_arr[3 * pos + 2] - pos_arr[3 * actual_cell_value + 2]) * (pos_arr[3 * pos + 2] - pos_arr[3 * actual_cell_value + 2]);
                                r2 = rx + ry + rz;

                                if (r2 <= kappa * kappa * h * h)
                                {
                                    // cout << "neighbour founded (before push)" << endl;
                                    // cout << "pos : " << pos << endl;
                                    // cout << "actual_cell_value : " << actual_cell_value << endl;
                                    neighbours_matrix[pos].push_back(actual_cell_value);
                                    // cout << "after first push_back" << endl;
                                    // cout << "len(neighbour_matrix) = " << neighbours_matrix.size() << endl;
                                    neighbours_matrix[actual_cell_value].push_back(pos);
                                    // cout << "neighbour founded (after push)" << endl;
                                }
                            }
                        }
                    }
                }
            }
        }

        cell_matrix[i_cell + j_cell * Nx + k_cell * Nx * Ny].erase(cell_matrix[i_cell + j_cell * Nx + k_cell * Nx * Ny].begin());
    }
}

void naiveAlgo(size_t nb_moving_part, vector<double> &pos_arr, vector<vector<int>> &neighbours_matrix, const double &h, const int &kappa)
{
    // Find neighbours for each particle
    for (int i = 0; i < nb_moving_part; i++)
    {
        for (int j = i + 1; j < nb_moving_part; j++)
        {
            double rx, ry, rz, r2;
            rx = (pos_arr[3 * i] - pos_arr[3 * j]) * (pos_arr[3 * i] - pos_arr[3 * j]);
            ry = (pos_arr[3 * i + 1] - pos_arr[3 * j + 1]) * (pos_arr[3 * i + 1] - pos_arr[3 * j + 1]);
            rz = (pos_arr[3 * i + 2] - pos_arr[3 * j + 2]) * (pos_arr[3 * i + 2] - pos_arr[3 * j + 2]);
            r2 = rx + ry + rz;
            if (r2 <= kappa * kappa * h * h)
            {
                neighbours_matrix[i].push_back(j);
                neighbours_matrix[j].push_back(i);
            }
        }
    }
}

void printNeighbours(vector<vector<int>> &neighbours_matrix_1, vector<vector<int>> &neighbours_matrix_2)
{
    for (int i = 0; i < neighbours_matrix_1.size(); i++)
    {
        std::cout << "Particle " << i << " : ";

        std::cout << "{";
        for (int j = 0; j < neighbours_matrix_1[i].size(); j++)
        {
            std::cout << neighbours_matrix_1[i][j];
            if (j != neighbours_matrix_1[i].size() - 1)
            {
                std::cout << ", ";
            }
        }
        std::cout << "} (Linked-list) VS {";

        for (int j = 0; j < neighbours_matrix_2[i].size(); j++)
        {
            std::cout << neighbours_matrix_2[i][j];
            if (j != neighbours_matrix_2[i].size() - 1)
            {
                std::cout << ", ";
            }
        }
        std::cout << "} (naive)\n \n";
    }
}
