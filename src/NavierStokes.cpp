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

    if (simParams.is_surface_tension){

        // Iterations over each particle
        #pragma omp parallel for
        for (int n = 0; n < nb_tot_part; n++){

            int size_neighbours = nb_neighbours[n];

            // Iterations over each associated neighbours 
            for (int idx = 0; idx < size_neighbours; idx++){

                int i_neig = neighbours[100*n + idx];

                vector<double> d_pos(3);
                for (int coord = 0; coord < 3; coord++)
                    d_pos[coord] = pos[3*n + coord] - pos[3*i_neig + coord];
                

                double r_ab = dist(pos, n, i_neig);
                double deriv = derive_cubic_spline(r_ab, geomParams, simParams);
                W[n][idx] = f_cubic_spline(r_ab, geomParams, simParams);;

                for (int coord = 0; coord < 3; coord++)
                    gradW[n][3*idx + coord] = (d_pos[coord] / r_ab) * deriv;
                
            }
        }
    }
    else{

        // Iterations over each particle
        #pragma omp parallel for
        for (int n = 0; n < nb_tot_part; n++){

            int size_neighbours = nb_neighbours[n];

            // Iterations over each associated neighbours 
            for (int idx = 0; idx < size_neighbours; idx++){

                int i_neig = neighbours[100*n + idx];

                vector<double> d_pos(3);
                for (int coord = 0; coord < 3; coord++)
                    d_pos[coord] = pos[3*n + coord] - pos[3*i_neig + coord];
                
                double r_ab = dist(pos, n, i_neig);
                double deriv = derive_wendland_quintic(r_ab, geomParams, simParams);
                W[n][idx] = f_wendland_quintic(r_ab, geomParams, simParams);;

                for (int coord = 0; coord < 3; coord++)
                    gradW[n][3*idx + coord] = (d_pos[coord] / r_ab) * deriv;
                
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
    int nb_part = simParams.nb_tot_part;

    #pragma omp parallel for
    for (int n = 0; n < nb_part; n++){

        if (state_equation == "Ideal gaz law") c[n] = c_0;        
        else if (state_equation == "Quasi incompresible fluid") c[n] = c_0 * pow(rho[n] / rho_0, 0.5 * (gamma - 1));
        else {
            cout << "Error : no state equation chosen" << endl;
            exit(1);
        } 

        if (c[n] < 0){
            cout << "Error, c ="<< c[n]<< " at timestep : " <<simParams.t << endl;
            exit(1);
        }
    
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

    #pragma omp parallel for
    for (int n = 0; n < nb_moving_part; n++){

        if (state_equation == "Ideal gaz law") p[n] =  (R * T / M) * (rho[n] / rho_0 - 1);

        else if (state_equation == "Quasi incompresible fluid"){
            double B = c_0 * c_0 * rho_0 / gamma;
            p[n] = B * (pow(rho[n] / rho_0, gamma) - 1);
        }
        else {
            cout << "Error : no state equation chosen" << endl;
            exit(1);
        }
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
        for (int n = 0; n < nb_moving_part; n++){

            int size_neighbours = nb_neighbours[n];

            // Iteration over each associated neighbours
            for (int idx = 0; idx < size_neighbours; idx++){

                vector<double> d_pos(3), d_u(3);
                int i_neig = neighbours[100*n + idx];

                for (int coord = 0; coord < 3; coord++){
                    d_pos[coord] = (pos[3*n + coord] - pos[3*i_neig + coord]);
                    d_u[coord] = (u[3*n + coord] - u[3*i_neig + coord]);
                }

                double c_a = c[n];
                double c_b = c[i_neig];
                double rho_a = rho[n];
                double rho_b = rho[i_neig];
                double c_ab = 0.5 * (c_a + c_b);
                double rho_ab = 0.5 * (rho_a + rho_b);
                double nu_2 = 0.01 * h * h;
                
                double u_ab_x_ab = dotProduct(d_u, d_pos);
                double x_ab_x_ab = dotProduct(d_pos, d_pos);
                double mu_ab = (h * u_ab_x_ab) / (x_ab_x_ab + nu_2);

                viscosity[n][idx] = (u_ab_x_ab < 0) ? 
                (-alpha * c_ab * mu_ab + beta * mu_ab * mu_ab) / rho_ab : 0;
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
    for (int n = 0; n < nb_part; n++){

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
                dot_product += (u_a - u_b) * gradW[n][3*idx + coord];
            }
            
            drhodt[n] += m_b * dot_product;
     
      
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
                      vector<double> &track_particle){


    int nb_moving_part = simParams.nb_moving_part;
    simParams.F_st_max = 0;

    // Compute pressure for all particles
    setPressure(geomParams, thermoParams, simParams, p, rho); 

    // Compute speed of sound for all particles
    setSpeedOfSound(geomParams, thermoParams, simParams, c, rho);

    // Compute artificial viscosity Π_ab for all particles
    setArtificialViscosity(geomParams, thermoParams, simParams, viscosity, 
                           neighbours, nb_neighbours, c, pos, rho, u); 

    vector<double> F_vol(3*simParams.nb_moving_part,0.0);
    
    
    if (simParams.is_surface_tension){
        InterfaceTrackingMath(simParams, geomParams, thermoParams,
                              nb_neighbours,neighbours, gradW,
                              mass,rho,type,pos, track_particle);

        surfaceTensionImprove(simParams, geomParams,thermoParams, nb_neighbours, neighbours, 
                              gradW, W, mass, rho, pos, F_vol, type, track_particle);
    }
    
    // Iterate over each particle
    #pragma omp parallel for
    for (int n = 0; n < nb_moving_part; n++){
        
        double rho_a = rho[n];
        double p_a = p[n];
        int size_neighbours = nb_neighbours[n];

        // Summation over neighbours
        for (int idx = 0; idx < size_neighbours; idx++){

            int i_neig = neighbours[100*n + idx];
            double pi_ab = viscosity[n][idx];
            double rho_b = rho[i_neig];
            double m_b = mass[i_neig];
            double p_b = p[i_neig];

            for (int coord = 0; coord < 3; coord++)
                dudt[3*n + coord] -= m_b * (p_b / (rho_b * rho_b) +
                                    p_a / (rho_a * rho_a) + pi_ab)* gradW[n][3*idx + coord];
            
            if (simParams.is_adhesion){
                
                double beta_ad = simParams.beta_adh;
                double r_ab = 0;
                vector<double> d_xyz(3);

                for (int coord = 0; coord < 3; coord++){
                
                    d_xyz[coord] = pos[3*n + coord] - pos[3*i_neig + coord];
                    r_ab += d_xyz[coord]*d_xyz[coord];
                }

                r_ab = sqrt(r_ab);
                double W_ab = W_adh(r_ab, geomParams, simParams);


                for (int coord = 0; coord < 3; coord++){
                    double boundary = 1.0 - type[i_neig];
                    F_vol[3*n + coord] -= beta_ad*boundary*mass[n]*m_b*W_ab*d_xyz[coord]/r_ab;
                }
            } 
        }
            
        double g = (simParams.is_gravity ? -9.81 : 0.0);
        double F_res = 0;

        for (int coord = 0; coord < 3; coord++){
            if(coord == 2)
                F_vol[3 * n + coord] += g;
            
            dudt[3 * n + coord] += F_vol[3 * n + coord];
            F_res += F_vol[3*n + coord]*F_vol[3*n + coord];
        }

        F_res = sqrt(F_res);
        simParams.F_st_max = (simParams.F_st_max > F_res ? simParams.F_st_max : F_res );     
    }

    if (simParams.PRINT) cout << "momentumEquation passed" << endl;
}

