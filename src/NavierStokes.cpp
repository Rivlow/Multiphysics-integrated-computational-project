#include <stdio.h>
#include <vector>
#include <omp.h>

#include "NavierStokes.h"
#include "find_neighbours.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"
#include "surface_tension.h"
using namespace std;


void computeGradW(GeomData &geomParams,    
                  SimulationData &simParams, 
                  vector<vector<double>> &gradW,
                  vector<vector<double>> &W,
                  vector<int> neighbours,
                  vector<double> nb_neighbours,
                  vector<double> pos){

    int nb_tot_part = simParams.nb_tot_part;

        // Iterations over each particle
        #pragma omp parallel for
        for (int i = 0; i < nb_tot_part; i++){

            int size_neighbours = nb_neighbours[i];

            // Iterations over each associated neighbours 
            for (int idx = 0; idx < size_neighbours; idx++){

                int j = neighbours[100*i + idx];

                vector<double> d_pos(3);
                for (int coord = 0; coord < 3; coord++)
                    d_pos[coord] = pos[3*i + coord] - pos[3*j + coord];
                

                double r_ij = dist(pos, i, j);
                double deriv = deriveCubicSpline(r_ij, geomParams, simParams);
                W[i][idx] = CubicSpline(r_ij, geomParams, simParams);
                //double deriv = deriveWendlandQuintic(r_ij, geomParams, simParams);
                //W[i][idx] = WendlandQuintic(r_ij, geomParams, simParams);                

                for (int coord = 0; coord < 3; coord++)
                    gradW[i][3*idx + coord] = (d_pos[coord] / r_ij) * deriv;
                
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
    int nb_part = simParams.nb_tot_part;


    if (state_equation == "Quasi incompresible fluid") {
        #pragma omp parallel for
        for (int i = 0; i < nb_part; i++){
             c[i] = c_0 * pow(rho[i] / rho_0, 0.5 * (gamma - 1));
        }
    }

    else if (state_equation == "Ideal gaz law"){; 
        #pragma omp parallel for
        for (int i = 0; i < nb_part; i++){
             c[i] = c_0;
        }
    }

    else {
        cout << "Error : no state equation chosen" << endl;
        exit(1);
    } 
    
        if (simParams.PRINT) cout << "setSpeedOfSound passed" << endl;
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
    int nb_moving_part = simParams.nb_moving_part;
    string state_equation = simParams.state_equation;


    if (state_equation == "Quasi incompresible fluid"){
        double B = c_0 * c_0 * rho_0 / gamma;
        #pragma omp parallel for
        for (int i = 0; i < nb_moving_part; i++){
            p[i] = B * (pow(rho[i] / rho_0, gamma) - 1);
        }
    }

   else if (state_equation == "Ideal gaz law"){
        #pragma omp parallel for
        for (int i = 0; i < nb_moving_part; i++){

         p[i] =  (R * T / M) * (rho[i] / rho_0 - 1);
        }

    }
    
    else {
        cout << "Error : no state equation chosen" << endl;
        exit(1);
    }

     if (simParams.PRINT) cout << "setPressure passed" << endl;
}


void setArtificialViscosity(GeomData &geomParams,    
                            ThermoData &thermoParams,
                            SimulationData &simParams, 
                            vector<vector<double>> &viscosity,
                            vector<int> &neighbours,
                            vector<double> &nb_neighbours,
                            vector<double> &c,
                            vector<double> &pos,
                            vector<double> &rho,
                            vector<double> &u){

    double beta = simParams.beta;
    double alpha = simParams.alpha;
    double h = geomParams.h;
    int nb_moving_part = simParams.nb_moving_part;
    int t = simParams.t;

    if (t > 0){

        // Iterations over each particle
        #pragma omp parallel for
        for (int i = 0; i < nb_moving_part; i++){

            int size_neighbours = nb_neighbours[i];

            // Iteration over each associated neighbours
            for (int idx = 0; idx < size_neighbours; idx++){

                vector<double> d_pos(3), d_u(3);
                int j = neighbours[100*i + idx];

                for (int coord = 0; coord < 3; coord++){
                    d_pos[coord] = (pos[3*i + coord] - pos[3*j + coord]);
                    d_u[coord] = (u[3*i + coord] - u[3*j + coord]);
                }

                double c_ij = 0.5 * (c[i] + c[j]);
                double rho_ij = 0.5 * (rho[i] + rho[j]);
                double nu_2 = 0.01 * h * h;
                
                double u_ij_x_ij = dotProduct(d_u, d_pos);
                double x_ij_x_ij = dotProduct(d_pos, d_pos);
                double mu_ab = (h * u_ij_x_ij) / (x_ij_x_ij + nu_2);

                viscosity[i][idx] = (u_ij_x_ij < 0) ? 
                (-alpha * c_ij * mu_ab + beta * mu_ab * mu_ab) / rho_ij : 0;
                
            }
        }
    }
    
    if (simParams.PRINT) cout << "setArtificialViscosity passed" << endl;
    
}

void continuityEquation(SimulationData &simParams, 
                        vector<int> &neighbours,
                        vector<double> &nb_neighbours,
                        vector<vector<double>> &gradW,
                        vector<double> &pos,
                        vector<double> &u,
                        vector<double> &drhodt,
                        vector<double> &rho,
                        vector<double> &mass){

    int nb_part = simParams.nb_tot_part;
             
    // Iterations over each particle
    #pragma omp parallel for
    for (int i = 0; i < nb_part; i++){

        int size_neighbours = nb_neighbours[i];

        // Summation over b = 1 -> nb_neighbours
        for (int idx = 0; idx < size_neighbours; idx++){

            int j = neighbours[100*i + idx];
            double dot_product = 0;
            
            // Dot product of u_ab with grad_a(W_ab)
            for (int coord = 0; coord < 3; coord++){

                double u_i = u[3*i + coord];
                double u_j = u[3*j + coord];
                dot_product += (u_i - u_j) * gradW[i][3*idx + coord];
            }
            
            drhodt[i] += mass[j] * dot_product;
        }
    }

    if (simParams.PRINT) cout << "continuityEquation passed" << endl;
}

void momentumEquation(GeomData &geomParams,    
                      ThermoData &thermoParams,
                      SimulationData &simParams, 
                      vector<int> &neighbours,
                      vector<double> &nb_neighbours,
                      vector<vector<double>> &gradW,
                      vector<vector<double>> W,
                      vector<vector<double>> &viscosity,
                      vector<double> &mass,
                      vector<double> &dudt,
                      vector<double> &rho,
                      vector<double> &p, 
                      vector<double> &c,
                      vector<double> &pos,
                      vector<double> &u,
                      vector<double> type,
                      vector<double> &colour,
                      vector<double> &R,
                      vector<double> &L,
                      vector<double> &i,
                      vector<double> &normal,
                      vector<double> &acc_vol,
                      vector<double> &track_particle,
                      vector<double> &Kappa,
                      vector<double> &dot_product){


    int nb_moving_part = simParams.nb_moving_part;
    simParams.acc_st_max = 0;

    // Compute pressure for all particles
    setPressure(geomParams, thermoParams, simParams, p, rho); 

    // Compute speed of sound for all particles
    setSpeedOfSound(geomParams, thermoParams, simParams, c, rho);

    // Compute artificial viscosity Π_ab for all particles
    setArtificialViscosity(geomParams, thermoParams, simParams, viscosity, 
                           neighbours, nb_neighbours, c, pos, rho, u); 
   
    if(simParams.is_surface_tension_1){
       
        surfaceTension( simParams, geomParams, thermoParams, nb_neighbours, neighbours, gradW, W,
                             mass, rho, pos, acc_vol, type, normal);
    }
    else if (simParams.is_surface_tension_2){
        
        InterfaceTrackingMath(simParams, geomParams, thermoParams,
                              nb_neighbours, neighbours, gradW,
                              mass, rho, type,pos, track_particle);

        surfaceTensionImprove(simParams, geomParams,thermoParams, nb_neighbours, neighbours, 
                              gradW, W, mass, rho, pos, type, colour, R, L, i, normal, acc_vol, track_particle, Kappa, dot_product);
    }

    // Iterate over each particle
    #pragma omp parallel for
    for (int i = 0; i < nb_moving_part; i++){
        
        double rho_i = rho[i];
        double p_i = p[i];
        int size_neighbours = nb_neighbours[i];

        // Summation over neighbours
        for (int idx = 0; idx < size_neighbours; idx++){

            int j = neighbours[100*i + idx];
            double pi_ij = viscosity[i][idx];
            double rho_j = rho[j];
            double m_j = mass[j];
            double p_j = p[j];
            
            for (int coord = 0; coord < 3; coord++){
                dudt[3*i + coord] -= m_j * (p_j / (rho_j * rho_j) +
                                    p_i / (rho_i * rho_i) + pi_ij)* gradW[i][3*idx + coord];
            }
                
            if (simParams.is_adhesion){
                
                double beta_ad = simParams.beta_adh;
                double r_ij = 0;
                vector<double> d_xyz(3);

                for (int coord = 0; coord < 3; coord++){
                
                    d_xyz[coord] = pos[3*i + coord] - pos[3*j + coord];
                    r_ij += d_xyz[coord]*d_xyz[coord];
                }

                r_ij = sqrt(r_ij);
                double W_ij = WAdh(r_ij, geomParams, simParams);

                double boundary = 1.0 - type[j];
                for (int coord = 0; coord < 3; coord++){
                    
                    acc_vol[3*i + coord] -= beta_ad*boundary*mass[i]*m_j*W_ij*d_xyz[coord]/r_ij;
                    
                }
            } 
        }
            
        double g = (simParams.is_gravity ? -9.81 : 0.0);
        double acc_res = 0;
        acc_vol[3*i + 2] += g;

        for (int coord = 0; coord < 3; coord++){
            
            dudt[3*i + coord] += acc_vol[3*i + coord];
            acc_res += acc_vol[3*i + coord]*acc_vol[3*i + coord];
        }
        
        acc_res = sqrt(acc_res);
        simParams.acc_st_max = (simParams.acc_st_max > acc_res ? simParams.acc_st_max : acc_res );     
    }

    if (simParams.PRINT) cout << "momentumEquation passed" << endl;
}

