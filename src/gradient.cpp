#include <stdio.h>
#include <vector>
#include "sorted_list.h"
#include "Kernel_functions.h"
using namespace std;



void gradW(vector<double> &pos_arr, vector<vector<unsigned>> &neighbours_matrix, vector<vector<double>> &gradW_matrix, 
                    double &h, unsigned &Nx, unsigned &Ny, unsigned &Nz){

    // Iterations over each particle
    for (unsigned pos = 0; pos < pos_arr.size()/3; pos++){

        vector<unsigned> &neighbours_list = neighbours_matrix[pos];
        vector<double> gradW_vect;

        // Iterations over each associated neighbours of prescribed particles
        for (unsigned idx_neighbour : neighbours_list){     
            
            double rx, ry, rz, r_ab;

            rx = (pos_arr[3*pos+0] - pos_arr[3*idx_neighbour+0])*(pos_arr[3*pos+0] - pos_arr[3*idx_neighbour+0]);
            ry = (pos_arr[3*pos+1] - pos_arr[3*idx_neighbour+1])*(pos_arr[3*pos+1] - pos_arr[3*idx_neighbour+1]);
            rz = (pos_arr[3*pos+2] - pos_arr[3*idx_neighbour+2])*(pos_arr[3*pos+2] - pos_arr[3*idx_neighbour+2]);
            r_ab = sqrt(rx + ry + rz);

            gradW_vect.push_back((pos_arr[3*pos+0] - pos_arr[3*idx_neighbour+0])/r_ab * derive_cubic_spline(r_ab, h));
            gradW_vect.push_back((pos_arr[3*pos+1] - pos_arr[3*idx_neighbour+1])/r_ab * derive_cubic_spline(r_ab, h));
            gradW_vect.push_back((pos_arr[3*pos+2] - pos_arr[3*idx_neighbour+2])/r_ab * derive_cubic_spline(r_ab, h));
        }

        gradW_matrix.push_back(gradW_vect);
        gradW_vect.clear();
    }
}

void setArtificialViscosity(vector<vector<double>> &artificial_visc,  vector<double> &pos_arr,  vector<vector<unsigned>> &neighbours_matrix, vector<double> &u_arr, 
                             double &c_ab,  double &rho_ab,  double &alpha,  double &beta,  double &h){


    vector<double> rel_displ(3), rel_velo(3);

    // Iterations over each particle
    for (size_t pos = 0; pos < pos_arr.size(); pos++){

        vector<unsigned> &neighbours_list = neighbours_matrix[pos];
        for (size_t idx_neighbour = 0; idx_neighbour < neighbours_list.size(); idx_neighbour++){   

            rel_displ[0] = (pos_arr[3*pos+0] - pos_arr[3*idx_neighbour+0]);
            rel_displ[1] = (pos_arr[3*pos+1] - pos_arr[3*idx_neighbour+1]);
            rel_displ[2] = (pos_arr[3*pos+2] - pos_arr[3*idx_neighbour+2]);

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



void stateEquation(double &p, double &c,  double &rho,  double &rho_0,  double &c_0,  double &R,  double &T,
                      double &M,  double &gamma,  string &state_equation_chosen){

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

void continuityEquation(vector<double> &pos_arr,  vector<double> &u_arr,  vector<vector<unsigned>> &neighbours_matrix, 
                         vector<vector<double>> &gradW_matrix, vector<double> &drhodt_arr, vector<double> &rho_arr,  vector<double> &mass_arr,  double &h){

    // Iterations over each particle                    
    for (size_t pos = 0; pos < pos_arr.size()/3; pos++){

         vector<unsigned> &neighbours_list = neighbours_matrix[pos];
         vector<double> &gradW_list = gradW_matrix[pos];

        double drhodt = 0, rho = 0;

        // Summation over b = 1 -> nb_neighbours
        for (size_t idx_neighbour = 0; idx_neighbour < neighbours_list.size(); idx_neighbour++){   

            // Dot product of u_ab with grad_a(W_ab)
            double dot_product = 0;
            for (size_t x = 0; x < 3; x++){

                dot_product += (u_arr[3*pos+x] - u_arr[3*neighbours_list[idx_neighbour]+x])*(gradW_list[idx_neighbour+x]);
            }

            double rx, ry, rz, r_ab;
            rx = (pos_arr[3*pos+0] - pos_arr[3*idx_neighbour+0])*(pos_arr[3*pos+0] - pos_arr[3*idx_neighbour+0]);
            ry = (pos_arr[3*pos+1] - pos_arr[3*idx_neighbour+1])*(pos_arr[3*pos+1] - pos_arr[3*idx_neighbour+1]);
            rz = (pos_arr[3*pos+2] - pos_arr[3*idx_neighbour+2])*(pos_arr[3*pos+2] - pos_arr[3*idx_neighbour+2]);
            r_ab = sqrt(rx + ry + rz);

            //rho += mass_arr[neighbours_list[idx_neighbour]]*f_cubic_spline(r_ab, h);
            drhodt += mass_arr[neighbours_list[idx_neighbour]]*dot_product;
        }

        //rho_arr[pos] = rho;
        drhodt_arr[pos] = drhodt;

       
    }
}


void momentumEquation(vector<vector<unsigned>> &neighbours_matrix,  vector<double> &mass_arr,  vector<vector<double>> &gradW_matrix, 
                      vector<double> &dudt_arr, vector<vector<double>> &artificial_visc,  vector<double> &rho_arr,  double &rho_0,  double &c_0,
                      vector<double> &p_arr,  double &R,  double &T,  double &M,  double &gamma, double &g, string &state_equation_chosen){

    // Iterations over each particle
    for (size_t pos = 0; pos < rho_arr.size(); pos++){
         vector<unsigned> &neighbours_list = neighbours_matrix[pos];
         vector<double> &gradW_list = gradW_matrix[pos];

        double p_a, c_a;
        stateEquation(p_a, c_a, rho_arr[pos], rho_0, c_0, R, T, M, gamma, state_equation_chosen);
        vector<double> dudt(3, 0.0), F_vol = {0.0, 0.0, mass_arr[pos]*g};

        // Summation over b = 1 -> nb_neighbours
        for (size_t idx_neighbour = 0; idx_neighbour < neighbours_list.size(); idx_neighbour++){   

            double p_b, c_b;
            stateEquation(p_b, c_b, rho_arr[neighbours_list[idx_neighbour]], rho_0, c_0, R, T, M, gamma, state_equation_chosen);
            double rho_a = rho_arr[pos], rho_b = rho_arr[neighbours_list[idx_neighbour]];
            double rho_ab = 0.5*(rho_a + rho_b);
            double pi_ab = 0 ; //[pos][neighbours_list[idx_neighbour]];
            double m_b = mass_arr[neighbours_list[idx_neighbour]];

            for (size_t cord = 0; cord < 3; cord++) {
                dudt[cord] += m_b*(p_b/(rho_b*rho_b) + p_a/(rho_a*rho_a) + pi_ab)*gradW_list[idx_neighbour+cord] + F_vol[cord];
            }
            
        }

        for (size_t cord = 0; cord < 3; cord++) {
            dudt[cord] *= -1; 
            dudt_arr[3*pos+cord] = dudt[cord];
        }

    }
}

