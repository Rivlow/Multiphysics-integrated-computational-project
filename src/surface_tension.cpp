#include <stdio.h>
#include <vector>
#include <cmath>
#include <omp.h>

#include "NavierStokes.h"
#include "find_neighbours.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"
#include "surface_tension.h"
#include "Kernel.h"

using namespace std;

/*
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
                    vector<double> &F_vol,
                    vector<double> type,
                    vector<double> normal_grad){

    vector<double> normal(3*simParams.nb_moving_part,0.0);
    //#pragma omp parallel for 
    for(int n = 0; n<simParams.nb_moving_part; n++){

        vector<double> &gradW = gradW_matrix[n];
        int size_neighbours = nb_neighbours[n];
        
        for(int idx = 0; idx<size_neighbours; idx++){ 

            int i_neig = neighbours[100*n + idx]; 
            double m_j = mass[i_neig];
            double rho_j = rho[i_neig];
            
            for( int coord = 0; coord <3; coord ++){

                double grad = gradW[3*idx+coord];
                normal[3*n+coord] += geomParams.h*m_j*grad/rho_j;
            }
        }
    }
    double alpha = simParams.alpha_st;

    //#pragma omp parallel for 
    for(int n = 0; n<simParams.nb_moving_part; n++){
    
        int size_neighbours = nb_neighbours[n];
        
        for(int idx = 0; idx<size_neighbours; idx++){
            
            if(type[neighbours[100*n + idx]] == 1){
                int i_neig = neighbours[100*n + idx];
                double K_ij = 2*thermoParams.rho_0/(rho[n]+rho[i_neig]);
                double r_ab = 0;
                vector<double> d_xyz(3);
                
                for (int coord = 0; coord < 3; coord++){
                    
                    d_xyz[coord] = pos[3 * n + coord] - pos[3 * i_neig + coord];
                    r_ab += d_xyz[coord]*d_xyz[coord];
                }
                
                r_ab = sqrt(r_ab);
                double W_ab = W_coh(r_ab,geomParams.kappa*geomParams.h, simParams);
                
                double m_a = mass[n];
                double m_b = mass[i_neig];
                double F_res = 0;
                
               // cout << "alpha*W_ab = " << alpha * W_ab << endl;
                for (int coord = 0; coord < 3; coord++){
                    
                    F_vol[3*n + coord] += -K_ij*(alpha * m_a * m_b * d_xyz[coord]*W_ab/r_ab 
                                    + alpha*(normal[3*n+coord]-normal[3*i_neig+coord]));
                
                    F_res += F_vol[3*n + coord]*F_vol[3*n + coord];
                }
                
                simParams.F_st_max = sqrt(F_res);
            }
        }     
    }  
}
*/


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
    for (int n = 0; n < simParams.nb_moving_part; n++){
    
        if (type[n] == 1){ // if moving particle

            int size_neighbours = nb_neighbours[n];
            double div_normal = 0;
            
            for (int idx = 0; idx < size_neighbours; idx++){

                int i_neig = neighbours[100*n + idx];
                if (type[i_neig] == 1){ // if moving neighbour

                    double dot_product = 0;
                    
                    for (int coord = 0; coord < 3; coord++)
                        dot_product += (pos[3*n + coord] - pos[3*i_neig + coord])*gradW[n][3*idx + coord];

                    double m_j = mass[i_neig];
                    double rho_j = rho[i_neig]; 
                    div_normal += (m_j/rho_j)*dot_product;
                }
            }  

            if (abs(div_normal) <= 1.7)
                track_particle[n] = 1;
        
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
                           vector<double> &F_vol,
                           vector<double> type,
                           vector<double> &track_particle){




    vector<int> N(simParams.nb_moving_part, 0);
    vector<double> colour(simParams.nb_moving_part, 0);
    vector<double> normal(3*simParams.nb_moving_part, 0);
    vector<double> R(simParams.nb_moving_part, 0);

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

    // Calculation of normal vector
    //#pragma omp parallel for
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

            if (track_particle[n] == 1 && track_particle[i_neig] == 0){  // case 1: part. i on free surface but not j

                for (int coord = 0; coord < 3; coord++)
                    normal[3*n+coord] -= (0 - colour[n])*(mass[n]/rho[n])*(-1*gradW[n][3*idx+coord]);
            }

            if (track_particle[n] == 0 && track_particle[i_neig] == 1){  // case 2 : part. j on free surface but not i

                //new gradient kernel has to be computed
                vector<double> d_pos(3);
                double r_ij = 0;

                for (int coord = 0; coord < 3; coord++){
                    
                    d_pos[coord] = 2*(pos[3*n + coord] - pos[3*i_neig + coord]); // vector ij' twice the lenght of classical ij vector
                    r_ij += d_pos[coord]*d_pos[coord];
                }

                r_ij = sqrt(r_ij);

                if (r_ij < geomParams.kappa*geomParams.h){ // if new imaginary j' particle in support domain of i

                    vector<double> new_gradW(3);
                    double deriv = derive_cubic_spline(r_ij, geomParams, simParams);
                
                    for (int coord = 0; coord < 3; coord++){
                        new_gradW[coord] = (d_pos[coord] / r_ij) * deriv;
                        normal[3*n+coord] -= (0 - colour[n])*(mass[n]/rho[n])*(new_gradW[coord]);

                    }
                }
            }            
        }


        // normalize normal 

        double norm = 0;

        for (int coord = 0; coord < 3; coord++)
            norm += normal[3*n+coord]*normal[3*n+coord];
        
        norm = sqrt(norm);
        R[n] = (norm <= 0.01/geomParams.h)? 0: 1;

        for (int coord = 0; coord < 3; coord++)
            normal[3*n+coord] /= (norm > 0)? norm : 1;
        
    }

  
    double Kappa = 0, L = 0;

    // Calculation of curvature
    //#pragma omp parallel for
    for (int n = 0; n < simParams.nb_moving_part; n++){

        int size_neighbours = nb_neighbours[n];

        for(int idx = 0; idx < size_neighbours; idx++){

            int i_neig = neighbours[100*n + idx];

            /*----------------------------*/
            /* classic normal calculation */
            /*----------------------------*/

            double dot_product = 0;
            for (int coord = 0; coord < 3; coord++)
                dot_product += (normal[3*i_neig + coord] - normal[3*n + coord])*gradW[n][3*idx + coord];
            
            Kappa += (R[n]*R[i_neig])*(mass[i_neig]/rho[i_neig])*dot_product;
            L += (R[n]*R[i_neig])*(mass[i_neig]/rho[i_neig])*W[n][idx];


            /*---------------------------------*/
            /* imaginary particle contribution */
            /*---------------------------------*/

            // case 1: part. i on free surface but not j
            if (track_particle[n] == 1 && track_particle[i_neig] == 0){ 

                vector<double> new_normal(3);
                double norm = 0;

                for (int coord = 0; coord < 3; coord++){
                    new_normal[coord] = 2*normal[3*n+coord] - normal[3*i_neig+coord]; 
                    norm += new_normal[coord]*new_normal[coord];
                }

                norm = sqrt(norm);
                double new_R = (norm <= 0.01/geomParams.h)? 0: 1;

                // normalize new vector
                for (int coord = 0; coord < 3; coord++)
                    new_normal[coord] = new_normal[coord]/norm;

                // contribution of normal ij'
                double dot_product = 0;

                for (int coord = 0; coord < 3; coord++)
                    dot_product += (new_normal[coord] - normal[3*n+coord])*(-1*gradW[n][3*idx+coord]);

                Kappa += (R[n]*new_R)*(mass[n]/rho[n])*dot_product;
                L += (R[n]*new_R)*(mass[n]/rho[n])*W[n][idx];
                
                continue;   
            }

            // case 2: part. j on free surface but not i
            if (track_particle[n] == 0 && track_particle[i_neig] == 1){ 

                //new normal and gradient kernel have to be computed
                vector<double> new_normal(3), d_pos(3);
                double norm = 0, r_ij = 0;
                
                for (int coord = 0; coord < 3; coord++){
                    new_normal[coord] = 2*normal[3*i_neig+coord] - normal[3*n+coord]; 
                    norm += new_normal[coord]*new_normal[coord];
                    d_pos[coord] = 2*(pos[3*n + coord] - pos[3*i_neig + coord]); // vector ij' twice the lenght of classical ij vector
                    r_ij += d_pos[coord]*d_pos[coord];
                } 

                norm = sqrt(norm);
                double new_R = (norm <= 0.01/geomParams.h)? 0: 1;

                // normalize new vector
                for (int coord = 0; coord < 3; coord++)
                    new_normal[coord] = new_normal[coord]/norm;

                // contribution of normal ij'
                double dot_product = 0;
                r_ij = sqrt(r_ij);

                if (r_ij < geomParams.kappa*geomParams.h){ // if new imaginary j' particle in support domain of i

                    vector<double> new_gradW(3);
                    double deriv = derive_cubic_spline(r_ij, geomParams, simParams);
                
                    for (int coord = 0; coord < 3; coord++){
                        new_gradW[coord] = (d_pos[coord] / r_ij) * deriv;
                        dot_product += (new_normal[coord] - normal[3*n+coord])*new_gradW[coord];
                    }
                }

                double new_W = f_cubic_spline(r_ij, geomParams, simParams);
                Kappa += (R[n]*new_R)*(mass[n]/rho[n])*dot_product;
                L += (R[n]*new_R)*(mass[n]/rho[n])*new_W;

            }            
        }

        Kappa *= -1;
        Kappa /= (L > 0)? L : 1; // correction of the truncated kernel support domain

        double F_res = 0;

        for (int coord = 0; coord < 3; coord++){

            F_vol[3*n+coord] += (thermoParams.sigma/rho[n])*Kappa*normal[3*n+coord];
            F_res += F_vol[3*n + coord]*F_vol[3*n + coord];
        }

        simParams.F_st_max = sqrt(F_res);

    }

    if (simParams.PRINT) cout << "surfaceTensionImprove passed" << endl;


}



