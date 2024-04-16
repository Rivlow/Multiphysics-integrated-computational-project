#include <stdio.h>
#include <vector>
#include <omp.h>

#include "gradient.h"
#include "find_neighbours.h"
#include "Kernel_functions.h"
#include "tools.h"
#include "structure.h"
#include "surface_tension.h"


using namespace std;

void gradW(SimulationData& params, 
           vector<vector<double>> &gradW_matrix,
           vector<vector<int>> &neighbours_matrix,
           vector<double> &pos){

    double h = params.h;
    int size_pos = pos.size()/3;
    int nb_moving_part = params.nb_moving_part;

    // Iterations over each particle
    #pragma omp parallel for
    for (int n = 0; n < nb_moving_part; n++){

        vector<int> &neighbours = neighbours_matrix[n];
        vector<double> &gradW = gradW_matrix[n];
        int size_neighbours = neighbours.size();

        //cout << "n : " << n << endl;
        // Iterations over each associated neighbours of prescribed particles
        for (int idx = 0; idx < size_neighbours; idx++){

            int i_neig = neighbours[idx];
            double r_ab, r_val = 0;
            vector<double> pos_val(3);

            for (int coord = 0; coord < 3; coord++){
            
                pos_val[coord] = pos[3 * n + coord] - pos[3 * i_neig + coord];
                r_val += pos_val[coord]*pos_val[coord];
            }

            r_ab = sqrt(r_val);
            double deriv = derive_cubic_spline(r_ab, h);
            /*cout << "idx : " << idx << endl;
            cout << "going to insert into gradW" << endl;*/
            
            for (int coord = 0; coord < 3; coord++){
                
                gradW[3 * idx + coord] = pos_val[coord] / r_ab * deriv;
            }
        }
    }

    if (params.PRINT){
            cout << "gradW passed" << endl;
    }
}

void setSpeedOfSound( SimulationData& params,
                     vector<double> &c,
                     vector<double> &rho){

    string state_equation = params.state_equation;
    double c_0 = params.c_0;
    double rho_0 = params.rho_0;
    double gamma = params.gamma;
    int size_rho = rho.size();

    #pragma omp parallel for
    for (int n = 0; n < size_rho; n++){

        if (state_equation == "Ideal gaz law"){
            c[n] = c_0;
        }

        if (state_equation == "Quasi incompresible fluid"){
            c[n] = c_0 * pow(rho[n] / rho_0, 0.5 * (gamma - 1));
        }
    }

}

void setPressure( SimulationData& params,
                 vector<double> &p,
                 vector<double> &rho){

    string state_equation = params.state_equation;
    double c_0 = params.c_0;
    double rho_0 = params.rho_0;
    double gamma = params.gamma;
    double R = params.R;
    double T = params.T;
    double M = params.M;
    bool PRINT = params.PRINT;
    int nb_moving_part = params.nb_moving_part;

    #pragma omp parallel for
    for (int n = 0; n < nb_moving_part; n++)
    {

        if (state_equation == "Ideal gaz law")
        {
            p[n] = (rho[n] / rho_0 - 1) * (R * T) / M;
        }

        if (state_equation == "Quasi incompresible fluid")
        {
            double B = c_0 * c_0 * rho_0 / gamma;
            p[n] = B * (pow(rho[n] / rho_0, gamma) - 1);
        }
    }

    if (PRINT){
            cout << "setPressure passed" << endl;
    }
}

void setArtificialViscosity( SimulationData& params,
                            int t,
                            vector<vector<double>> &artificial_visc_matrix,
                            vector<vector<int>> &neighbours_matrix,
                            vector<double> &c,
                            vector<double> &pos,
                            vector<double> &rho,
                            vector<double> &u){

    double beta = params.beta;
    double alpha = params.alpha;
    double h = params.h;
    bool PRINT = params.PRINT;
    int nb_moving_part = params.nb_moving_part;


    if (t == 0){
        #pragma omp parallel for
        for (int n = 0; n < nb_moving_part; n++){

            vector<int> &neighbours = neighbours_matrix[n];
            int size_neighbours = neighbours.size();

            for (int idx_neighbour = 0; idx_neighbour < size_neighbours; idx_neighbour++){
                artificial_visc_matrix[n].push_back(0.0);
            }
        }
    }

    else{

        vector<double> rel_displ(3), rel_vel(3);
        int size_pos = pos.size()/3;


        // Iterations over each particle
        #pragma omp parallel for
        for (int n = 0; n < size_pos; n++){

            vector<int> &neighbours = neighbours_matrix[n];
            int size_neighbours = neighbours.size();

            // Iteration over each associated neighbours
            for (int idx = 0; idx < size_neighbours; idx++){

                int i_neig = neighbours[idx];

                rel_displ[0] = (pos[3 * n + 0] - pos[3 * i_neig + 0]);
                rel_displ[1] = (pos[3 * n + 1] - pos[3 * i_neig + 1]);
                rel_displ[2] = (pos[3 * n + 2] - pos[3 * i_neig + 2]);

                rel_vel[0] = (u[3 * n + 0] - u[3 * i_neig + 0]);
                rel_vel[1] = (u[3 * n + 1] - u[3 * i_neig + 1]);
                rel_vel[2] = (u[3 * n + 2] - u[3 * i_neig + 2]);

                double u_ab_x_ab = 0, x_ab_2 = 0;

                // Dot product
                for (int cord = 0; cord < 3; cord++){
                    u_ab_x_ab += rel_vel[cord] * rel_displ[cord];
                    x_ab_2 += rel_displ[cord] * rel_displ[cord];
                }

                double c_a = c[n];
                double c_b = c[i_neig];
                double rho_a = rho[n];
                double rho_b = rho[i_neig];

                double c_ab = 0.5 * (c_a + c_b);
                double rho_ab = 0.5 * (rho_a + rho_b);
                double nu_2 = 0.01 * h * h;
                double mu_ab = (h * u_ab_x_ab) / (x_ab_2 + nu_2);

                artificial_visc_matrix[n].push_back((u_ab_x_ab < 0) ? (-alpha * c_ab * mu_ab + beta * mu_ab * mu_ab) / rho_ab : 0);
            }
        }
    }
    
    if (PRINT){
            cout << "setArtificialViscosity passed" << endl;
    }
}

void continuityEquation( SimulationData& params,
                        vector<vector<int>> &neighbours_matrix,
                        vector<vector<double>> &gradW_matrix,
                        vector<double> &pos,
                        vector<double> &u,
                        vector<double> &drhodt,
                        vector<double> &rho,
                        vector<double> &mass){

    bool PRINT = params.PRINT;
    int size_pos = pos.size()/3;
             
    // Iterations over each particle
    #pragma omp parallel for
    for (int n = 0; n < size_pos; n++){

        vector<int> &neighbours = neighbours_matrix[n];
        vector<double> &gradW= gradW_matrix[n];
        int size_neighbours =neighbours.size();

        // Summation over b = 1 -> nb_neighbours
        for (int idx = 0; idx < size_neighbours; idx++){

            
            double dot_product = 0;
            int i_neig = neighbours[idx];
            double m_b = mass[i_neig];

            // Dot product of u_ab with grad_a(W_ab)
            for (int cord = 0; cord < 3; cord++){

                double u_a = u[3 * n + cord];
                double u_b = u[3 * i_neig + cord];
                dot_product += (u_a - u_b) * gradW[3*idx + cord];
            }
            
            drhodt[n] += m_b * dot_product;
        }

    }

    if (PRINT){
            cout << "continuityEquation passed" << endl;
    }
}

void momentumEquation( SimulationData& params,
                      int t,
                      vector<vector<int>> &neighbours_matrix,
                      vector<vector<double>> &gradW_matrix,
                      vector<vector<double>> &artificial_visc_matrix,
                      vector<double> &mass,
                      vector<double> &dudt,
                      vector<double> &rho,
                      vector<double> &p, 
                      vector<double> &c,
                      vector<double> &pos,
                      vector<double> &u){


    double g = params.g;
    bool PRINT = params.PRINT;
    int nb_moving_part = params.nb_moving_part;
    vector<double> F_vol(pos.size());

    // Compute pressure for all particles
    setPressure(params, p, rho); 

    // Compute speed of sound for all particles
    setSpeedOfSound(params, c, rho);

    // Compute artificial viscosity Π_ab for all particles
    setArtificialViscosity(params, t, artificial_visc_matrix, neighbours_matrix,
                            c, pos, rho, u); 

    // Iterations over each particle
    //printArray(F_vol,F_vol.size(), "F_vol");
    //surfaceTension( params, gradW_matrix, neighbours_matrix, mass, rho,
                    //pos, F_vol);
    //printArray(F_vol,F_vol.size(), "F_vol");

    #pragma omp parallel for
    for (int n = 0; n < nb_moving_part; n++){

        vector<int> &neighbours = neighbours_matrix[n];
        vector<double> &gradW = gradW_matrix[n];
        vector<double> &artificial_visc = artificial_visc_matrix[n];
        
        double rho_a = rho[n];
        double p_a = p[n];

        for (int cord = 0; cord < 3; cord++){

            // Summation over b = 1 -> nb_neighbours
            for (int idx = 0; idx < int(neighbours.size()); idx++){

                int i_neig = neighbours[idx];
                double pi_ab = artificial_visc[idx];
                double rho_b = rho[i_neig];
                double m_b = mass[i_neig];
                double p_b = p[i_neig];

                dudt[3 * n + cord] += m_b * (p_b / (rho_b * rho_b) + p_a / (rho_a * rho_a) + pi_ab) 
                                            * gradW[3*idx + cord];
            }

            dudt[3 * n + cord] *= -1;
            dudt[3 * n + cord] += F_vol[3*n+cord];
        }
    }
    
    if (PRINT){
        cout << "momentumEquation passed" << endl;
    }
}


