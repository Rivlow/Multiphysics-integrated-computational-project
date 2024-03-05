#include <stdio.h>
#include <vector>
#include "sorted_list.h"
#include "Kernel_functions.h"
using namespace std;

void gradW(const vector<double> &part_pos, const vector<vector<unsigned>> &neighbours_matrix, vector<vector<double>> &gradW_matrix, 
           double L[3], const double &h, const int &Nx, const int &Ny, const int &Nz){

    // Iterations over each particle
    for (unsigned pos = 0; pos < part_pos.size()/3; pos++){

        const vector<unsigned> &neighbours_list = neighbours_matrix[pos];
        cout << "nb neighbours : " << neighbours_list.size() << " \n" << endl;
        vector<double> gradW_vect;

        // Iterations over each associated neighbours of prescribed particles
        for (unsigned idx_neighbour : neighbours_list){     
            
            double rx, ry, rz, r_ab;

            rx = (part_pos[3*pos+0] - part_pos[3*idx_neighbour+0])*(part_pos[3*pos+0] - part_pos[3*idx_neighbour+0]);
            ry = (part_pos[3*pos+1] - part_pos[3*idx_neighbour+1])*(part_pos[3*pos+1] - part_pos[3*idx_neighbour+1]);
            rz = (part_pos[3*pos+2] - part_pos[3*idx_neighbour+2])*(part_pos[3*pos+2] - part_pos[3*idx_neighbour+2]);
            r_ab = sqrt(rx + ry + rz);

            gradW_vect.push_back((part_pos[3*pos+0] - part_pos[3*idx_neighbour+0])/r_ab * derive_cubic_spline(r_ab, h));
            gradW_vect.push_back((part_pos[3*pos+1] - part_pos[3*idx_neighbour+1])/r_ab * derive_cubic_spline(r_ab, h));
            gradW_vect.push_back((part_pos[3*pos+2] - part_pos[3*idx_neighbour+2])/r_ab * derive_cubic_spline(r_ab, h));
        }

        gradW_matrix.push_back(gradW_vect);
        gradW_vect.clear();
    }

    /*
    for (size_t i = 0; i < gradW_matrix.size(); i++) {
            std::cout << "Pour le " << i << "e element : (";
            
            for (size_t j = 0; j < gradW_matrix[i].size(); j++) {
                std::cout << gradW_matrix[i][j];
                
                if (j < gradW_matrix[i].size() - 1) {
                    std::cout << ", ";
                }
            }

            std::cout << ")" << std::endl;
        }
    */
    
}

void continuityEquation(const vector<double> &part_pos, const vector<vector<unsigned>> &neighbours_matrix, 
                        const vector<vector<double>> &gradW_matrix, vector<double> &drhodt_arr, vector<double> &rho_arr, const double &mass, const double &h){

    for (size_t pos = 0; pos < part_pos.size()/3; pos++){

        const vector<unsigned> &neighbours_list = neighbours_matrix[pos];
        const vector<double> &gradW_list = gradW_matrix[pos];

        double drhodt = 0, rho = 0;

        // Summation over b = 1 -> nb_neighbours
        for (size_t idx_neighbour = 0; idx_neighbour < neighbours_list.size(); idx_neighbour++){   

            /* !!!!!! A CHECK : est ce que "part_pos[3*neighbours_list[idx_neighbour]+x]" et "gradW_list[idx_neighbour+x]" 
            font bien reference au meme au meme voisin ???? Normalement oui mais pas sur a 100% "*/

            // Dot product of u_ab with grad_a(W_ab)


            /* ATTENTION !!!! Ici j'ai mis "(x,y,z)" plutot que "(u_x,u_y,u_z)"" uniquement pour voir ce que ca affichait -> a changer plus tard*/
            double dot_product = 0;
            for (size_t x = 0; x < 3; x++){

                dot_product += (part_pos[3*pos+x] - part_pos[3*neighbours_list[idx_neighbour]+x])*(gradW_list[idx_neighbour+x]);
            }

            double rx, ry, rz, r_ab;
            rx = (part_pos[3*pos+0] - part_pos[3*idx_neighbour+0])*(part_pos[3*pos+0] - part_pos[3*idx_neighbour+0]);
            ry = (part_pos[3*pos+1] - part_pos[3*idx_neighbour+1])*(part_pos[3*pos+1] - part_pos[3*idx_neighbour+1]);
            rz = (part_pos[3*pos+2] - part_pos[3*idx_neighbour+2])*(part_pos[3*pos+2] - part_pos[3*idx_neighbour+2]);
            r_ab = sqrt(rx + ry + rz);

            rho += mass*f_cubic_spline(r_ab, h);
            drhodt += mass*dot_product; // mass of particles is constant ?? Or rather use mass_matrix (which would be set initially) in case of ?

        }

        rho_arr[pos] = rho;
        drhodt_arr[pos] = drhodt;
    }
}

void momentumEquation(const vector<vector<unsigned>> &neighbours_matrix, const double &mass, const vector<vector<double>> &gradW_matrix, const vector<double> &rho_arr, const double &rho_0, const double &c_0,
                      vector<double> &p_arr, const double &R, const double &T, const double &M, const double &gamma, const string &state_equation_chosen){

    for (size_t a = 0; a < rho_arr.size(); a++){

        double p = stateEquation(rho_arr[a], rho_0, c_0, R, T, M, gamma, state_equation_chosen);

       

            const vector<unsigned> &neighbours_list = neighbours_matrix[a];
            const vector<double> &gradW_list = gradW_matrix[a];

            double dudt = 0;

            // Summation over b = 1 -> nb_neighbours
            for (size_t idx_neighbour = 0; idx_neighbour < neighbours_list.size(); idx_neighbour++){   

         

   
            }

        .
    }
}


double stateEquation(const double &rho, const double &rho_0, const double &c_0, const double &R, const double &T,
                     const double &M, const double &gamma, const string state_equation_chosen){

    if (state_equation_chosen == "Ideal gaz law"){
        double p = (rho/rho_0 - 1)*(rho*R*T)/M;
    }
    if (state_equation_chosen == "Quasi incompresible fluid"){
        double B = c_0*c_0*rho_0/gamma;
        double p = B*(pow(rho/rho_0, gamma) - 1);
    }
}