#include <stdio.h>
#include <vector>
#include <omp.h>

#include "gradient.h"
#include "find_neighbours.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"
#include "surface_tension.h"


using namespace std;

void gradW(GeomData &geomParams,    
           SimulationData &simParams, 
           vector<vector<double>> &gradW_matrix,
           vector<vector<double>> &W_matrix,
           vector<int> &neighbours,
           vector<double> &nb_neighbours,
           vector<double> &pos){


    double h = geomParams.h;
    int nb_part = simParams.nb_part;

    // Iterations over each particle
    #pragma omp parallel for
    for (int n = 0; n < nb_part; n++){

        vector<double> &gradW = gradW_matrix[n];
        int size_neighbours = nb_neighbours[n];

        // Iterations over each associated neighbours 
        for (int idx = 0; idx < size_neighbours; idx++){

            int i_neig = neighbours[100*n + idx];
            double r_ab = 0;
            vector<double> d_xyz(3);

            for (int coord = 0; coord < 3; coord++){
                
                d_xyz[coord] = pos[3 * n + coord] - pos[3 * i_neig + coord];
                r_ab += d_xyz[coord]*d_xyz[coord];
            }

            r_ab = sqrt(r_ab);
            double deriv = derive_cubic_spline(r_ab, h);
            double W = f_cubic_spline(r_ab, h);
            W_matrix[n][idx] = W;


            for (int coord = 0; coord < 3; coord++){
                gradW[3 * idx + coord] = d_xyz[coord] / r_ab * deriv;
            }
        }
    }

    if (simParams.PRINT)cout << "gradW passed" << endl;
    
}


void setSpeedOfSound(GeomData &geomParams,    
                     ThermoData &thermoParams,
                     SimulationData &simParams, 
                     vector<double> &c,
                     vector<double> &rho){

    string state_equation = simParams.state_equation;
    double c_0 = thermoParams.c_0;
    double rho_0 = thermoParams.rho_0;
    double gamma = thermoParams.gamma;
    bool PRINT = simParams.PRINT;
    int nb_part = simParams.nb_part;

    #pragma omp parallel for
    for (int n = 0; n < nb_part; n++){

        if (state_equation == "Ideal gaz law") c[n] = c_0;        
        else if (state_equation == "Quasi incompresible fluid") c[n] = c_0 * pow(rho[n] / rho_0, 0.5 * (gamma - 1));
        else {
            cout << "Error : no state equation chosen" << endl;
            exit(1);
        } 

        if (c[n] < 0){
            cout << "c ="<< c[n]<< " at timestep : " <<simParams.t << endl;
            exit(1);
        }
    
    }

        if (PRINT) cout << "setSpeedOfSound passed" << endl;
}

void setPressure(GeomData &geomParams,    
                 ThermoData &thermoParams,
                 SimulationData &simParams, 
                 vector<double> &p,
                 vector<double> &rho){

    double c_0 = thermoParams.c_0;
    double rho_0 = thermoParams.rho_0;
    double gamma = thermoParams.gamma;
    double R = thermoParams.R;
    double T = thermoParams.T;
    double M = thermoParams.M;
    bool PRINT = simParams.PRINT;
    int nb_moving_part = simParams.nb_moving_part;
    string state_equation = simParams.state_equation;

    #pragma omp parallel for
    for (int n = 0; n < nb_moving_part; n++){

        if (state_equation == "Ideal gaz law") p[n] = (rho[n] / rho_0 - 1) * (R * T) / M;

        else if (state_equation == "Quasi incompresible fluid"){
            double B = c_0 * c_0 * rho_0 / gamma;
            p[n] = B * (pow(rho[n] / rho_0, gamma) - 1);
        }
        else {
            cout << "Error : no state equation chosen" << endl;
            exit(1);
        }
    }

    if (PRINT) cout << "setPressure passed" << endl;
    
}

void setArtificialViscosity(GeomData &geomParams,    
                            ThermoData &thermoParams,
                            SimulationData &simParams, 
                            vector<vector<double>> &pi_matrix,
                            vector<int> &neighbours,
                            vector<double> &nb_neighbours,
                            vector<double> &c,
                            vector<double> &pos,
                            vector<double> &rho,
                            vector<double> &u){

    double beta = simParams.beta;
    double alpha = simParams.alpha;
    double h = geomParams.h;
    bool PRINT = simParams.PRINT;
    int nb_moving_part = simParams.nb_moving_part;
    int t = simParams.t;

    if (t == 0){
        #pragma omp parallel for
        for (int n = 0; n < nb_moving_part; n++){

            int size_neighbours = nb_neighbours[n];

            for (int idx = 0; idx < size_neighbours; idx++)
                pi_matrix[n][idx] = 0;
        }
    }

    else{

        vector<double> rel_displ(3), rel_vel(3);

        // Iterations over each particle
        #pragma omp parallel for
        for (int n = 0; n < nb_moving_part; n++){

            int size_neighbours = nb_neighbours[n];

            // Iteration over each associated neighbours
            for (int idx = 0; idx < size_neighbours; idx++){

                int i_neig = neighbours[100*n + idx];

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

                pi_matrix[n][idx] = (u_ab_x_ab < 0) ? 
                (-alpha * c_ab * mu_ab + beta * mu_ab * mu_ab) / rho_ab : 0;
            }
        }
    }
    
    if (PRINT) cout << "setArtificialViscosity passed" << endl;
    
}

void continuityEquation(SimulationData& simParams,
                        vector<int> &neighbours,
                        vector<double> &nb_neighbours,
                        vector<vector<double>> &gradW_matrix,
                        vector<double> &pos,
                        vector<double> &u,
                        vector<double> &drhodt,
                        vector<double> &rho,
                        vector<double> &mass){

    bool PRINT = simParams.PRINT;
    int nb_part = simParams.nb_part;
             
    // Iterations over each particle
    #pragma omp parallel for
    for (int n = 0; n < nb_part; n++){

        vector<double> &gradW = gradW_matrix[n];
        int size_neighbours = nb_neighbours[n];

        // Summation over b = 1 -> nb_neighbours
        for (int idx = 0; idx < size_neighbours; idx++){

            double dot_product = 0;
            int i_neig = neighbours[100*n + idx];
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

    if (PRINT) cout << "continuityEquation passed" << endl;
}

void momentumEquation(GeomData &geomParams,    
                      ThermoData &thermoParams,
                      SimulationData &simParams, 
                      vector<int> &neighbours,
                      vector<double> &nb_neighbours,
                      vector<int> &track_surface,
                      vector<double> &N_smoothed,
                      vector<vector<double>> &gradW_matrix,
                      vector<vector<double>> W_matrix,
                      vector<vector<double>> &pi_matrix,
                      vector<double> &mass,
                      vector<double> &dudt,
                      vector<double> &rho,
                      vector<double> &p, 
                      vector<double> &c,
                      vector<double> &pos,
                      vector<double> &u){



    double g = (simParams.is_gravity)? -9.81 : 0;
    bool PRINT = simParams.PRINT;
    int nb_moving_part = simParams.nb_moving_part;

    // Compute pressure for all particles
    setPressure(geomParams, thermoParams, simParams, p, rho); 

    // Compute speed of sound for all particles
    setSpeedOfSound(geomParams, thermoParams, simParams, c, rho);

    // Compute artificial viscosity Π_ab for all particles
    setArtificialViscosity(geomParams, thermoParams, simParams, pi_matrix, 
                           neighbours, nb_neighbours, c, pos, rho, u); 

    vector<double> F_vol(3*simParams.nb_moving_part,0.0);

    if (simParams.is_surface_tension)
        surfaceTension(simParams, geomParams,thermoParams, nb_neighbours, neighbours, 
                       track_surface, N_smoothed, gradW_matrix, W_matrix, mass, rho, pos, F_vol);



    // Iterate over each particle
    #pragma omp parallel for
    for (int n = 0; n < nb_moving_part; n++){

        vector<double> &gradW = gradW_matrix[n];
        vector<double> &artificial_visc = pi_matrix[n];
        
        double rho_a = rho[n];
        double p_a = p[n];

        for (int cord = 0; cord < 3; cord++){

            int size_neighbours = nb_neighbours[n];

            // Summation over b = 1 -> nb_neighbours
            for (int idx = 0; idx < size_neighbours; idx++){

                int i_neig = neighbours[100*n + idx];
                double pi_ab = artificial_visc[idx];
                double rho_b = rho[i_neig];
                double m_b = mass[i_neig];
                double p_b = p[i_neig];

                dudt[3 * n + cord] += m_b * (p_b / (rho_b * rho_b) +
                                      p_a / (rho_a * rho_a) + pi_ab)* gradW[3*idx + cord];
            }

            dudt[3 * n + cord] *= -1;
            dudt[3 * n + cord] += F_vol[3 * n + cord];

            if(cord == 2)
                dudt[3 * n + cord] += g;
            
        }
    }

    if (PRINT) cout << "momentumEquation passed" << endl;
}


