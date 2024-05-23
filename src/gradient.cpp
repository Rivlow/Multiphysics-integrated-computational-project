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
           vector<int> neighbours,
           vector<double> nb_neighbours,
           vector<double> pos,
           vector<double> rho,
           vector<double> &normal,
           vector<double> mass){


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
            double deriv = derive_cubic_spline(r_ab, h, simParams);
            double W = f_cubic_spline(r_ab, h, simParams);
            double m_b = mass[i_neig];
            double rho_b = rho[i_neig];
            W_matrix[n][idx] = W;
            double h = geomParams.h;

            for (int coord = 0; coord < 3; coord++){
                gradW[3 * idx + coord] = (d_xyz[coord] / r_ab) * deriv;
                normal[3 * n + coord] += gradW[3 * idx + coord]*h*m_b/rho_b;
            }
        }
    }

    if (simParams.PRINT)cout << "gradW passed" << endl;
    //printArray(normal, normal.size(),"normal  1");
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
                for (int coord = 0; coord < 3; coord++){
                    u_ab_x_ab += rel_vel[coord] * rel_displ[coord];
                    x_ab_2 += rel_displ[coord] * rel_displ[coord];
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
            for (int coord = 0; coord < 3; coord++){

                double u_a = u[3 * n + coord];
                double u_b = u[3 * i_neig + coord];
                dot_product += (u_a - u_b) * gradW[3*idx + coord];
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
                      vector<vector<double>> &gradW_matrix,
                      vector<vector<double>> W_matrix,
                      vector<vector<double>> &pi_matrix,
                      vector<double> &mass,
                      vector<double> &dudt,
                      vector<double> &rho,
                      vector<double> &p, 
                      vector<double> &c,
                      vector<double> &pos,
                      vector<double> &u,
                      vector<double> type,
                      vector<double> normal){


    double g = (simParams.is_gravity ? -9.81 : 0.0);
    bool PRINT = simParams.PRINT;
    int nb_moving_part = simParams.nb_moving_part;
    simParams.F_st_max = 0;
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
                        gradW_matrix, W_matrix, mass, rho, pos, F_vol,type, normal);

    //printArray(F_vol,F_vol.size(),"fvol0");
    double alpha = simParams.alpha_st;
    // Iterate over each particle
    #pragma omp parallel for
    for (int n = 0; n < nb_moving_part; n++){

        vector<double> &gradW = gradW_matrix[n];
        vector<double> &artificial_visc = pi_matrix[n];
        
        double rho_a = rho[n];
        double p_a = p[n];
        int size_neighbours = nb_neighbours[n];

            // Summation over b = 1 -> nb_neighbours
        for (int idx = 0; idx < size_neighbours; idx++){

            int i_neig = neighbours[100*n + idx];
            double pi_ab = artificial_visc[idx];
            double rho_b = rho[i_neig];
            double m_b = mass[i_neig];
            double p_b = p[i_neig];
            for (int coord = 0; coord < 3; coord++){
                dudt[3 * n + coord] += m_b * (p_b / (rho_b * rho_b) +
                                    p_a / (rho_a * rho_a) + pi_ab)* gradW[3*idx + coord];
            }
            
            /*if (simParams.is_surface_tension) {
                    
                double K_ij = 2*thermoParams.rho_0/(rho[n]+rho[i_neig]);
                double r_ab = 0;
                vector<double> d_xyz(3);
            
                for (int coord = 0; coord < 3; coord++){
                    
                    d_xyz[coord] = pos[3 * n + coord] - pos[3 * i_neig + coord];
                    r_ab += d_xyz[coord]*d_xyz[coord];
                }
                double h = geomParams.h;
                double kappa = geomParams.kappa;
                r_ab = sqrt(r_ab);
                double W_ab = W_coh(r_ab,kappa*h, simParams);
                //cout << W_ab << endl;
                double m_a = mass[n];
                double m_b = mass[i_neig];
                double boundary =  type[i_neig];
                
                for (int coord = 0; coord < 3; coord++){
                    
                    F_vol[3*n + coord] += -boundary*K_ij*(alpha * m_a * m_b * d_xyz[coord]*W_ab/r_ab 
                                    + alpha*(normal[3*n+coord]-normal[3*i_neig+coord])); 
                    

                }
            }*/

            if(simParams.is_adhesion){
                
                double beta_ad = simParams.beta_adh;
                double r_ab = 0;
                vector<double> d_xyz(3);

                for (int coord = 0; coord < 3; coord++){
                
                    d_xyz[coord] = pos[3 * n + coord] - pos[3 * i_neig + coord];
                    r_ab += d_xyz[coord]*d_xyz[coord];
                }

                r_ab = sqrt(r_ab);
                double kh = geomParams.kappa*geomParams.h;
                double W_ab = W_adh(r_ab, kh, simParams);


                for (int coord = 0; coord < 3; coord++){
                    double boundary = 1.0 - type[i_neig];
                    F_vol[3 * n + coord] -= beta_ad*boundary*mass[n]*m_b*W_ab*d_xyz[coord]/r_ab;
                }
            } 
        }
            
        double F_res = 0;
        for (int coord = 0; coord < 3; coord++){
            if(coord == 2)
                F_vol[3 * n + coord] += g;
            
            dudt[3 * n + coord] *= -1;
            dudt[3 * n + coord] += F_vol[3 * n + coord];
            F_res += F_vol[3 * n + coord]*F_vol[3 * n + coord];
        }
        F_res = sqrt(F_res);
        simParams.F_st_max = (simParams.F_st_max > F_res ? simParams.F_st_max : F_res );
            
    }
    

    //printArray(F_vol, F_vol.size(), "F_vol");
    if (PRINT) cout << "momentumEquation passed" << endl;
}


