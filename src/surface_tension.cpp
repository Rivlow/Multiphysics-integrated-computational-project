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
    for(int i = 0; i<simParams.nb_moving_part; i++){

        vector<double> &gradW = gradW_matrix[i];
        int size_neighbours = nb_neighbours[i];
        
        for(int idx = 0; idx<size_neighbours; idx++){ 

            int j = neighbours[100*i + idx]; 
            double m_j = mass[j];
            double rho_j = rho[j];
            
            for( int coord = 0; coord <3; coord ++){

                double grad = gradW[3*idx+coord];
                normal[3*i+coord] += geomParams.h*m_j*grad/rho_j;
            }
        }
    }
    double alpha = simParams.alpha_st;

    //#pragma omp parallel for 
    for(int i = 0; i<simParams.nb_moving_part; i++){
    
        int size_neighbours = nb_neighbours[i];
        
        for(int idx = 0; idx<size_neighbours; idx++){
            
            if(type[neighbours[100*i + idx]] == 1){
                int j = neighbours[100*i + idx];
                double K_ij = 2*thermoParams.rho_0/(rho[i]+rho[j]);
                double r_ab = 0;
                vector<double> d_xyz(3);
                
                for (int coord = 0; coord < 3; coord++){
                    
                    d_xyz[coord] = pos[3 * i + coord] - pos[3 * j + coord];
                    r_ab += d_xyz[coord]*d_xyz[coord];
                }
                
                r_ab = sqrt(r_ab);
                double W_ab = W_coh(r_ab,geomParams.kappa*geomParams.h, simParams);
                
                double m_a = mass[i];
                double m_b = mass[j];
                double F_res = 0;
                
               // cout << "alpha*W_ab = " << alpha * W_ab << endl;
                for (int coord = 0; coord < 3; coord++){
                    
                    F_vol[3*i + coord] += -K_ij*(alpha * m_a * m_b * d_xyz[coord]*W_ab/r_ab 
                                    + alpha*(normal[3*i+coord]-normal[3*j+coord]));
                
                    F_res += F_vol[3*i + coord]*F_vol[3*i + coord];
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
                           vector<double> &F_vol,
                           vector<double> type,
                           vector<double> &colour,
                           vector<double> &R,
                           vector<double> &N,
                           vector<double> &normal,
                           vector<double> &track_particle,
                           vector<double> &Kappa,
                           vector<double> &dot_product){


    // Calculation of i 
    #pragma omp parallel for 
    for (int i = 0; i < simParams.nb_moving_part; i++){

        int size_neighbours = nb_neighbours[i];

        for(int idx = 0; idx < size_neighbours; idx++){

            int j = neighbours[100*i + idx];

            if (track_particle[j]){ // if boundary part. in neighbourhood -> assigned 1 to part. "i"
                N[i] = 1;
                continue;
            }
        }
    }


    // Calculation of colour function
    #pragma omp parallel for
    for (int i = 0; i < simParams.nb_moving_part; i++){

        int size_neighbours = nb_neighbours[i];
        if (N[i]){ // free surface particule in support domain

            for(int idx = 0; idx < size_neighbours; idx++){

                int j = neighbours[100*i + idx];
                colour[i] += (mass[j]/rho[j])*W[i][idx]; // c_j^0 = 1 for fluid 
                
            }
        }
            
        else 
            colour[i] = 1;
            
    }

    // Calculation of normal vector
    #pragma omp parallel for
    for (int i = 0; i < simParams.nb_moving_part; i++){

        int size_neighbours = nb_neighbours[i];

        for(int idx = 0; idx < size_neighbours; idx++){

            int j = neighbours[100*i + idx];

            /*----------------------------*/
            /* classic normal calculation */
            /*----------------------------*/

            for (int coord = 0; coord < 3; coord++)
                normal[3*i+coord] -= (colour[j] - colour[i]) * (mass[j] / rho[j]) * gradW[i][3*idx+coord];

            /*---------------------------------*/
            /* imaginary particle contribution */
            /*---------------------------------*/

            if (track_particle[i] == 1 && track_particle[j] == 0){  // case 1: part. i on free surface but not j

                for (int coord = 0; coord < 3; coord++)
                    normal[3*i+coord] -= (0 - colour[i])*(mass[i]/rho[i])*(-1*gradW[i][3*idx+coord]);
            }

            if (track_particle[i] == 0 && track_particle[j] == 1){  // case 2 : part. j on free surface but not i

                //new gradient kernel has to be computed
                vector<double> d_pos(3);
                double r_ij = 0;

                for (int coord = 0; coord < 3; coord++){
                    
                    d_pos[coord] = 2*(pos[3*i + coord] - pos[3*j + coord]); // vector ij' twice the lenght of classical ij vector
                    r_ij += d_pos[coord]*d_pos[coord];
                }

                r_ij = sqrt(r_ij);

                // if new imaginary j' particle in support domain of i
                if (r_ij < geomParams.kappa*geomParams.h){ 
                    vector<double> new_gradW(3);
                    double deriv = derive_cubic_spline(r_ij, geomParams, simParams);
                
                    for (int coord = 0; coord < 3; coord++){
                        new_gradW[coord] = (d_pos[coord] / r_ij) * deriv;
                        normal[3*i+coord] -= (0 - colour[i])*(mass[i]/rho[i])*(new_gradW[coord]);

                    }
                }
            }          
        }

        // normalize normal 
        double norm = 0;

        for (int coord = 0; coord < 3; coord++)
            norm += normal[3*i+coord]*normal[3*i+coord];
        
        norm = sqrt(norm);
        R[i] = (norm <= 0.01/geomParams.h)? 0: 1;

        for (int coord = 0; coord < 3; coord++)
            normal[3*i+coord] /= (norm > 0)? norm : 1;
        
    }

    // Calculation of curvature
    #pragma omp parallel for
    for (int i = 0; i < simParams.nb_moving_part; i++){

        double L = 0;
        int size_neighbours = nb_neighbours[i];

        for(int idx = 0; idx < size_neighbours; idx++){

            int j = neighbours[100*i + idx];

            /*----------------------------*/
            /* classic normal calculation */
            /*----------------------------*/
            dot_product[i] = 0;

            for (int coord = 0; coord < 3; coord++)
                dot_product[i] += (normal[3*j + coord] - normal[3*i + coord]) * gradW[i][3*idx + coord];

            Kappa[i] -= (R[i]*R[j]) * (mass[j]/rho[j]) * dot_product[i];
            L += (R[i]*R[j]) * (mass[j]/rho[j]) * W[i][idx];


            /*---------------------------------*/
            /* imaginary particle contribution */
            /*---------------------------------*/

            // case 1: part. i on free surface but not j
            if (track_particle[i] == 1 && track_particle[j] == 0){ 

                vector<double> new_normal(3);
                double norm = 0;

                for (int coord = 0; coord < 3; coord++){
                    new_normal[coord] = 2 * normal[3*i+coord] - normal[3*j+coord]; 
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
                    dot_product += (new_normal[coord] - normal[3*i+coord]) * (-1*gradW[i][3*idx+coord]);

                Kappa[i] -= (R[i]*new_R) * (mass[i]/rho[i]) * dot_product;
                L += (R[i]*new_R) * (mass[i]/rho[i]) * W[i][idx];
                
                continue;   
            }

            // case 2: part. j on free surface but not i
            if (track_particle[i] == 0 && track_particle[j] == 1){ 

                //new normal and  have to be computed
                vector<double> new_normal(3), d_pos(3);
                double norm = 0, r_ij = 0;
                
                for (int coord = 0; coord < 3; coord++){
                    new_normal[coord] = 2 * normal[3*j+coord] - normal[3*i+coord]; 
                    norm += new_normal[coord] * new_normal[coord];
                    d_pos[coord] = 2*(pos[3*i + coord] - pos[3*j + coord]); // vector ij' twice the lenght of classical ij vector
                    r_ij += d_pos[coord] * d_pos[coord];
                } 

                norm = sqrt(norm);
                double new_R = (norm <= 0.01/geomParams.h)? 0: 1;

                // normalize new vector n_ij'
                for (int coord = 0; coord < 3; coord++)
                    new_normal[coord] = new_normal[coord]/norm;

                r_ij = sqrt(r_ij);

                if (r_ij < geomParams.kappa*geomParams.h){ // if part. j'  in support domain of i

                    // new kernel and gradient kernel have to be computed
                    vector<double> new_gradW(3);
                    double dot_product = 0;
                    double deriv = derive_cubic_spline(r_ij, geomParams, simParams);
                    double new_W = f_cubic_spline(r_ij, geomParams, simParams);
                
                    for (int coord = 0; coord < 3; coord++){
                        new_gradW[coord] = (d_pos[coord] / r_ij) * deriv;
                        dot_product += (new_normal[coord] - normal[3*i+coord])*new_gradW[coord];
                    }

                    Kappa[i] -= (R[i]*new_R) * (mass[i]/rho[i]) * dot_product;
                    L += (R[i]*new_R) * (mass[i]/rho[i]) * new_W;

                }
            }            
        }

        Kappa[i] /= (L > 0)? L : 1; // correction of the truncated kernel support domain
        double F_res = 0;

        for (int coord = 0; coord < 3; coord++){

            F_vol[3*i+coord] += (thermoParams.sigma/rho[i])*Kappa[i]*normal[3*i+coord];
            F_res += F_vol[3*i + coord]*F_vol[3*i + coord];
        }

        simParams.F_st_max = sqrt(F_res);

    }

    if (simParams.PRINT) cout << "surfaceTensionImprove passed" << endl;


}



