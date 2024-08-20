#include <stdio.h>
#include <vector>
#include <cmath>
#include <omp.h>
#include "export.h"

#include "NavierStokes.h"
#include "find_neighbours.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"
#include "surface_tension.h"
#include "Kernel.h"

using namespace std;


void surfaceTension(SimulationData& simParams,
                    GeomData &geomParams,
                    ThermoData &thermoParams,
                    vector<double> nb_neighbours,
                    vector<int> neighbours,
                    vector<vector<double>> gradW_matrix,
                    vector<vector<double>> W_matrix,
                    vector<double> mass,
                    vector<double> rho,
                    vector<double> pos,
                    vector<double> &acc_vol,
                    vector<double> type,
                    vector<double> normal){

   
    // compute normal
    #pragma omp parallel for 
    for(int i = 0; i<simParams.nb_moving_part; i++){

        vector<double> &gradW = gradW_matrix[i];
        int size_neighbours = nb_neighbours[i];
        
        for(int idx = 0; idx<size_neighbours; idx++){ 

            int j = neighbours[100*i + idx]; 
            double m_j = mass[j];
            double rho_j = rho[j];

            for(int coord = 0; coord <3; coord ++){

                double grad = gradW[3*idx+coord];
                normal[3*i+coord] += geomParams.h*m_j*grad/rho_j;
            } 
        }
    }

    double alpha = simParams.alpha_st;

    // Compute surface tension forces
    #pragma omp parallel for 
    for(int i = 0; i<simParams.nb_moving_part; i++){
    
        int size_neighbours = nb_neighbours[i];
        
        for(int idx = 0; idx<size_neighbours; idx++){
            
            if(type[neighbours[100*i + idx]] == 1){

                int j = neighbours[100*i + idx];
                double K_ij = 2*thermoParams.rho_0/(rho[i]+rho[j]);
                double r_ij = 0;
                vector<double> d_xyz(3);
                
                for (int coord = 0; coord < 3; coord++){
                    
                    d_xyz[coord] = pos[3*i + coord] - pos[3*j + coord];
                    r_ij += d_xyz[coord]*d_xyz[coord];
                }
                
                r_ij = sqrt(r_ij);
                double W_ij = WCoh(r_ij, geomParams, simParams);
                double m_i = mass[i];
                double m_j = mass[j];
                double acc_res = 0;
                
                for (int coord = 0; coord < 3; coord++){
                    
                    acc_vol[3*i + coord] += -K_ij*(alpha * m_i*m_j * (d_xyz[coord]/r_ij)*W_ij 
                                    + alpha * m_i * (normal[3*i+coord]-normal[3*j+coord]));
                    acc_res += acc_vol[3*i + coord]*acc_vol[3*i + coord];
                }
                
                simParams.acc_st_max = sqrt(acc_res);
            }
        }     
    }  
}



void InterfaceTrackingMath(SimulationData simParams,
                           GeomData geomParams,
                           ThermoData thermoParams,
                           vector<double> nb_neighbours,
                           vector<int> neighbours,
                           vector<vector<double>> gradW,
                           vector<double> mass,
                           vector<double> rho,
                           vector<double> type,
                           vector<double> pos,
                           vector<double> &track_particle){
    


    
    #pragma omp parallel for 
    for (int i = 0; i < simParams.nb_moving_part; i++){
    
        if (type[i] == 1){ // if moving particle

            int size_neighbours = nb_neighbours[i];
            double div_normal = 0;
            
            for (int idx = 0; idx < size_neighbours; idx++){

                int j = neighbours[100*i + idx];
                if (type[j] == 1){ // if moving neighbour

                    double dot_product = 0;
                    
                    for (int coord = 0; coord < 3; coord++)
                        dot_product += (pos[3*i + coord] - pos[3*j + coord])*gradW[i][3*idx + coord];

                    double m_j = mass[j];
                    double rho_j = rho[j]; 
                    div_normal += (m_j/rho_j)*dot_product;
                }
            }  

            if (abs(div_normal) <= 1.7)
                track_particle[i] = 1;
        
        }  
    }

    if (simParams.PRINT) cout << "InterfaceTrackingMath passed" << endl;

}



void surfaceTensionImprove(SimulationData& simParams,
                           GeomData &geomParams,
                           ThermoData &thermoParams,
                           vector<double> &nb_neighbours,
                           vector<int> &neighbours,
                           vector<vector<double>> &gradW,
                           vector<vector<double>> &W,
                           vector<double> &mass,
                           vector<double> &rho,
                           vector<double> &pos,
                           vector<double> type,
                           vector<double> &colour,
                           vector<double> &R,
                           vector<double> &L,
                           vector<double> &N,
                           vector<double> &normal,
                           vector<double> &acc_vol,
                           vector<double> &track_particle,
                           vector<double> &Kappa,
                           vector<double> &dot_product){


    // Calculation of N 
    #pragma omp parallel for 
    for (int n = 0; n < simParams.nb_moving_part; n++){

        int size_neighbours = nb_neighbours[n];

        for(int idx = 0; idx < size_neighbours; idx++){

            int i_neig = neighbours[100*n + idx];

            if (track_particle[i_neig]){ // if boundary part. in neighbourhood -> assigned 1 to part. "n"
                N[n] = 1;
                continue;
            }
        }
    }


    // Calculation of colour function
    #pragma omp parallel for
    for (int n = 0; n < simParams.nb_moving_part; n++){

        int size_neighbours = nb_neighbours[n];
        if (N[n]){ // free surface particule in support domain

            for(int idx = 0; idx < size_neighbours; idx++){

                int i_neig = neighbours[100*n + idx];
                colour[n] += (mass[i_neig]/rho[i_neig])*W[n][idx]; // c_j^0 = 1 for fluid (to be checked)
                
            }
        }
            
        else 
            colour[n] = 1;
            
    }

    vector<double> normal_normali(3*simParams.nb_moving_part,0.0);

    // Calculation of normal vector
    vector<double> imaginary_part1(simParams.nb_moving_part,0.0);
    #pragma omp parallel for
    for (int n = 0; n < simParams.nb_moving_part; n++){

        int size_neighbours = nb_neighbours[n];

        for(int idx = 0; idx < size_neighbours; idx++){

            int i_neig = neighbours[100*n + idx];

            /*----------------------------*/
            /* classic normal calculation */
            /*----------------------------*/

            for (int coord = 0; coord < 3; coord++)
                normal[3*n+coord] -= (colour[i_neig] - colour[n]) * (mass[i_neig] / rho[i_neig]) * gradW[n][3*idx+coord];

            
            /*---------------------------------*/
            /* imaginary particle contribution */
            /*---------------------------------*/

            
            // case 1: part. i on free surface but not j
            if (track_particle[n] == 1 && track_particle[i_neig] == 0){  
                imaginary_part1[n] +=1;
                for (int coord = 0; coord < 3; coord++)
                    normal[3*n+coord] -= (0 - colour[n])*(mass[n]/rho[n])*(-1*gradW[n][3*idx+coord]);
            }


              // case 2 : part. j on free surface but not i
            if (track_particle[n] == 0 && track_particle[i_neig] == 1){
                
                //new gradient kernel has to be computed
                vector<double> d_pos(3);
                double r_ij = 0;

                for (int coord = 0; coord < 3; coord++){
                    
                    d_pos[coord] = 2*(pos[3*n + coord] - pos[3*i_neig + coord]); // vector ij' twice the lenght of classical ij vector
                    r_ij += d_pos[coord]*d_pos[coord];
                }

                r_ij = sqrt(r_ij);

                if (r_ij < geomParams.kappa*geomParams.h){ // if new imaginary j' particle in support domain of i
                    imaginary_part1[n] +=1;
                    vector<double> new_gradW(3);
                    double deriv = deriveWendlandQuintic(r_ij, geomParams, simParams);
                
                    for (int coord = 0; coord < 3; coord++){
                        new_gradW[coord] = (d_pos[coord] / r_ij) * deriv;
                        normal[3*n+coord] -= (0 - colour[n])*(mass[n]/rho[n])*(new_gradW[coord]);

                    }
                }
            }         
        }


        // Normalize normal 
        double norm = 0;

        for (int coord = 0; coord < 3; coord++)
            norm += normal[3*n+coord]*normal[3*n+coord];
        
        norm = sqrt(norm);
        R[n] = (norm <= 0.01/geomParams.h)? 0: 1;
        
        for (int coord = 0; coord < 3; coord++)
            normal_normali[3*n+coord] = (norm > 0)? normal[3*n+coord]/norm : normal[3*n+coord];
        
    }

    vector<double> imaginary_part2(simParams.nb_moving_part,0.0);

    // Calculation of curvature
    #pragma omp parallel for
    for (int n = 0; n < simParams.nb_moving_part; n++){

        double L = 0;
        int size_neighbours = nb_neighbours[n];

        for(int idx = 0; idx < size_neighbours; idx++){

            int i_neig = neighbours[100*n + idx];

            /*----------------------------*/
            /* classic normal calculation */
            /*----------------------------*/
            dot_product[n] = 0;

            for (int coord = 0; coord < 3; coord++)
                dot_product[n] += (normal_normali[3*i_neig + coord] - normal_normali[3*n + coord]) * gradW[n][3*idx + coord];

            Kappa[n] -= (R[n]*R[i_neig]) * (mass[i_neig]/rho[i_neig]) * dot_product[n];
            L += (R[n]*R[i_neig]) * (mass[i_neig]/rho[i_neig]) * W[n][idx];


            /*---------------------------------*/
            /* imaginary particle contribution */
            /*---------------------------------*/
            
            // case 1: part. i on free surface but not j
            if (track_particle[n] == 1 && track_particle[i_neig] == 0){ 

                imaginary_part2[n] +=1;
                vector<double> new_normal(3);
                double norm = 0;

                for (int coord = 0; coord < 3; coord++){
                    new_normal[coord] = 2 * normal_normali[3*n+coord] - normal_normali[3*i_neig+coord]; 
                    norm += new_normal[coord] * new_normal[coord];
                }

                norm = sqrt(norm);
                double new_R = (norm <= 0.01/geomParams.h)? 0: 1;

                // normalize new vector
                for (int coord = 0; coord < 3; coord++)
                    new_normal[coord] = new_normal[coord]/norm;

                // contribution of normal ij'
                double dot_product = 0;

                for (int coord = 0; coord < 3; coord++)
                    dot_product += (new_normal[coord] - normal_normali[3*n+coord]) * (-1*gradW[n][3*idx+coord]);

                Kappa[n] -= (R[n]*new_R) * (mass[n]/rho[n]) * dot_product;
                L += (R[n]*new_R) * (mass[n]/rho[n]) * W[n][idx];
                
                
                continue;   
            }

            // case 2: part. j on free surface but not i
            if (track_particle[n] == 0 && track_particle[i_neig] == 1){ 
                
                //new normal and gradient kernel have to be computed
                vector<double> new_normal(3), d_pos(3);
                double norm = 0, r_ij = 0;
                
                for (int coord = 0; coord < 3; coord++){
                    new_normal[coord] = 2 * normal_normali[3*i_neig+coord] - normal_normali[3*n+coord]; 
                   
                    norm += new_normal[coord] * new_normal[coord];
                    d_pos[coord] = 2*(pos[3*n + coord] - pos[3*i_neig + coord]); // vector ij' twice the lenght of classical ij vector
                    r_ij += d_pos[coord] * d_pos[coord];
                } 

                norm = sqrt(norm);
                double new_R = (norm <= 0.01/geomParams.h)? 0: 1;

                // normalize new vector
                for (int coord = 0; coord < 3; coord++)
                    new_normal[coord] = new_normal[coord]/norm;

                // contribution of normal ij'
                double dot_product = 0;
                r_ij = sqrt(r_ij);

                if (r_ij < geomParams.kappa*geomParams.h){ // if part. j'  in support domain of i
                    imaginary_part2[n] +=1;
                    vector<double> new_gradW(3);
                    double deriv = deriveWendlandQuintic(r_ij, geomParams, simParams);
                
                    for (int coord = 0; coord < 3; coord++){
                        new_gradW[coord] = (d_pos[coord] / r_ij) * deriv;
                        dot_product += (new_normal[coord] - normal_normali[3*n+coord])*new_gradW[coord];
                    }
                    double new_W = WendlandQuintic(r_ij, geomParams, simParams);
                    Kappa[n] -= (R[n]*new_R) * (mass[n]/rho[n]) * dot_product;
                    L += (R[n]*new_R) * (mass[n]/rho[n]) * new_W;
                }
            }   
        }

        Kappa[n] /= (L > 0)? L : 1; // correction of the truncated kernel support domain
        double acc_res = 0;

        for (int coord = 0; coord < 3; coord++){

            acc_vol[3*n+coord] += (thermoParams.sigma/rho[n])*Kappa[n]*normal_normali[3*n+coord];
            acc_res += acc_vol[3*n + coord]*acc_vol[3*n + coord];
        }

        simParams.acc_st_max = sqrt(acc_res);

    }
    
    if (simParams.PRINT) cout << "surfaceTensionImprove passed" << endl;

}
