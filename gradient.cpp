#include <stdio.h>
#include <vector>
#include "sorted_list.h"
#include "Kernel_functions.h"
using namespace std;

void gradW(const vector<double> &part_pos, const vector<vector<unsigned>> &neighbours_matrix, vector<vector<double>> &gradW_matrix, 
           double L[3], const double &h, const int &Nx, const int &Ny, const int &Nz){

    for (unsigned pos = 0; pos < part_pos.size()/3; pos++){

        const vector<unsigned> &neighbours_list = neighbours_matrix[pos];
        vector<double> gradW_vect;

        for (unsigned idx_neighbour : neighbours_list){     
            
            double rx, ry, rz, r_ab;

            cout << "pos : " << pos << " and idx_neighbour : " << idx_neighbour << endl;
            rx = (part_pos[3*pos] - part_pos[3*idx_neighbour])*(part_pos[3*pos] - part_pos[3*idx_neighbour]);
            ry = (part_pos[3*pos + 1] - part_pos[3*idx_neighbour+1])*(part_pos[3*pos + 1] - part_pos[3*idx_neighbour+1]);
            rz = (part_pos[3*pos + 2] - part_pos[3*idx_neighbour+2])*(part_pos[3*pos + 2] - part_pos[3*idx_neighbour+2]);
            r_ab = rx + ry + rz;

            cout << " relative x-distance : " << part_pos[3*pos] - part_pos[3*idx_neighbour]<< "\n ";
            cout << "relative y-distance : " << part_pos[3*pos+1] - part_pos[3*idx_neighbour+1]<< "\n ";
            cout << "relative z-distance : " << part_pos[3*pos+2] - part_pos[3*idx_neighbour+2] << "\n " << endl;

            cout << "r_ab/h = " << r_ab/h << "\n "<< endl;

            gradW_vect.push_back((part_pos[3*pos+0] - part_pos[3*idx_neighbour+0])/r_ab * derive_cubic_spline(r_ab, h));
            gradW_vect.push_back((part_pos[3*pos+1] - part_pos[3*idx_neighbour+1])/r_ab * derive_cubic_spline(r_ab, h));
            gradW_vect.push_back((part_pos[3*pos+2] - part_pos[3*idx_neighbour+2])/r_ab * derive_cubic_spline(r_ab, h));
        }

        gradW_matrix.push_back(gradW_vect);
    }

        
    for (size_t i = 0; i < gradW_matrix.size(); ++i) {
            std::cout << "Pour le " << i << "e element : (";
            
            // Parcours de chaque élément dans le vecteur interne
            for (size_t j = 0; j < gradW_matrix[i].size(); ++j) {
                std::cout << gradW_matrix[i][j];
                
                // Ajouter une virgule et un espace pour tous les éléments sauf le dernier
                if (j < gradW_matrix[i].size() - 1) {
                    std::cout << ", ";
                }
            }

            std::cout << ")" << std::endl;
        }
    
}

void continuityEquation(const vector<double> &part_pos, const vector<vector<unsigned>> &neighbours_matrix, 
                        const vector<vector<double>> &gradW_matrix, vector<double> &drhodt_arr, const double &mass){

    for (unsigned pos = 0; pos < part_pos.size()/3; pos++){

        const vector<unsigned> &neighbours_list = neighbours_matrix[pos];
        const vector<double> &gradW_list = gradW_matrix[pos];
        vector<double> gradW_vect;

        double drhodt = 0;

        // Summation over b = 1 -> N
        for (unsigned idx_neighbour : neighbours_list){     
            
            // Dot product of u_ab with grad_a(W_ab)
            double dot_product;
            for (unsigned x = 0; x < 3; x++){
                dot_product += (part_pos[3*pos+x] - part_pos[3*idx_neighbour+x])*(gradW_list[3*pos*x]);
            }

            drhodt += mass*dot_product; // mass of particles is constant ?? Or rather use mass_matrix (which would be set initially) in case of ?

        }

        drhodt_arr[pos] = drhodt;

        
    }
    
}