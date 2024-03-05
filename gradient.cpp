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

double setArtificialViscosity(vector<vector<double>> &artificial_visc, const vector<double> &part_pos, const vector<vector<unsigned>> &neighbours_matrix, vector<double> &u_arr, 
                            const double &c_ab, const double &rho_ab, const double &alpha, const double &beta, const double &h){


    vector<double> rel_displ(3), rel_velo(3);

    // Iterations over each particle
    for (size_t pos = 0; pos < part_pos.size(); pos++){

        const vector<unsigned> &neighbours_list = neighbours_matrix[pos];
        for (size_t idx_neighbour = 0; idx_neighbour < neighbours_list.size(); idx_neighbour++){   

            rel_displ[0] = (part_pos[3*pos+0] - part_pos[3*idx_neighbour+0]);
            rel_displ[1] = (part_pos[3*pos+1] - part_pos[3*idx_neighbour+1]);
            rel_displ[2] = (part_pos[3*pos+2] - part_pos[3*idx_neighbour+2]);

            rel_velo[0] = (u_arr[3*pos+0] - u_arr[3*idx_neighbour+0]);
            rel_velo[1] = (u_arr[3*pos+1] - u_arr[3*idx_neighbour+1]);
            rel_velo[2] = (u_arr[3*pos+2] - u_arr[3*idx_neighbour+2]);

            double res = 0, xa_xb = 0;
            // Dot product
            for (size_t idx = 0; idx < 3; idx++){
                res += rel_velo[idx]*rel_displ[idx];
                xa_xb += rel_displ[idx]*rel_displ[idx];
            }

            double nu_2 = 0.01*h*h;
            double mu_ab = (h*res)/(xa_xb + nu_2);
            artificial_visc[pos].push_back((res < 0) ? (-alpha*c_ab*mu_ab + beta*mu_ab*mu_ab)/rho_ab : 0);

        }
    }
}

void continuityEquation(const vector<double> &part_pos, const vector<vector<unsigned>> &neighbours_matrix, 
                        const vector<vector<double>> &gradW_matrix, vector<double> &drhodt_arr, vector<double> &rho_arr, const double &mass, const double &h){

    // Iterations over each particle                    
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
            drhodt += mass*dot_product;

        }

        rho_arr[pos] = rho;
        drhodt_arr[pos] = drhodt;
    }
}

void momentumEquation(const vector<vector<unsigned>> &neighbours_matrix, const double &mass, const vector<vector<double>> &gradW_matrix, const vector<double> &rho_arr, const double &rho_0, const double &c_0,
                      vector<double> &p_arr, const double &R, const double &T, const double &M, const double &gamma, const string &state_equation_chosen){

    // Iterations over each particle
    for (size_t pos = 0; pos < rho_arr.size(); pos++){

        const vector<unsigned> &neighbours_list = neighbours_matrix[pos];
        const vector<double> &gradW_list = gradW_matrix[pos];

        double p_a, c_a;
        stateEquation(p_a, c_a, rho_arr[pos], rho_0, c_0, R, T, M, gamma, state_equation_chosen);

        vector<double> dudt(3);

        // Summation over b = 1 -> nb_neighbours
        for (size_t idx_neighbour = 0; idx_neighbour < neighbours_list.size(); idx_neighbour++){   

            double p_b, c_b;
            stateEquation(p_b, c_b, rho_arr[neighbours_list[idx_neighbour]], rho_0, c_0, R, T, M, gamma, state_equation_chosen);

            double rho_a = rho_arr[pos], rho_b = rho_arr[neighbours_list[idx_neighbour]];
            double rho_ab = 0.5*(rho_a + rho_b);
            //double pi_ab = se

            dudt[0] = m_b*(p_b/(rho_b*rho_b) + p_a/(rho_a*rho_a) + pi_ab)*gradW_list[idx_neighbour];

            

        }
    }
}


void stateEquation(double &p, double &c, const double &rho, const double &rho_0, const double &c_0, const double &R, const double &T,
                     const double &M, const double &gamma, const string state_equation_chosen){

    if (state_equation_chosen == "Ideal gaz law"){
        p = (rho/rho_0 - 1)*(rho*R*T)/M;
        c = c_0;
    }
    if (state_equation_chosen == "Quasi incompresible fluid"){
        double B = c_0*c_0*rho_0/gamma;
        p = B*(pow(rho/rho_0, gamma) - 1);
        c = c_0*pow(rho/rho_0, 0.5*(gamma-1));
    }

}