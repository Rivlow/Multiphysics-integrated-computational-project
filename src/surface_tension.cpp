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
                        dot_product += (pos[3 * n + coord] - pos[3 * i_neig + coord])*gradW[n][3*idx+coord];

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
    vector<double> color(simParams.nb_moving_part, 0);
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

    //printArray(N, N.size(), "N");



    // Calculation of color function
    #pragma omp parallel for
    for (int n = 0; n < simParams.nb_moving_part; n++){

        int size_neighbours = nb_neighbours[n];
        if (N[n]){ // free surface particule in support domain

            
            for(int idx = 0; idx < size_neighbours; idx++){

                int i_neig = neighbours[100*n + idx];
                color[n] += (mass[i_neig]/rho[i_neig])*W[n][idx]; // c_j^0 = 1 for fluid (to be checked)
                
            }
        }
            
        else 
            color[n] = 1;
            
    }


    // Calculation of normal vector
    //#pragma omp parallel for
    for (int n = 0; n < simParams.nb_moving_part; n++){

        int size_neighbours = nb_neighbours[n];

        for(int idx = 0; idx < size_neighbours; idx++){

            int i_neig = neighbours[100*n + idx];

            // classic normal calculation 
            for (int coord = 0; coord < 3; coord++)
                normal[3*n+coord] -= (color[i_neig] - color[n]) * (mass[i_neig] / rho[i_neig]) * gradW[n][3*idx+coord];
            
            // imaginary particle contribution
            if (track_particle[n] == 1 && track_particle[i_neig] == 0){ // particle on free surface but not neighbour 


                for (int coord = 0; coord < 3; coord++){
                    normal[3*n+coord] -= (0 - color[n])*(mass[n]/rho[n])*(-1*gradW[n][3*idx+coord]);

                }
                                   
            }

            if (N[n]){ // at least one free surface neighbour in support domain
                if (track_particle[n] == 0 && track_particle[i_neig] == 1){ // some neighbours on free surface but not particle i 

                    // new gradient kernel has to be computed
                    vector<double> d_xyz(3);
                    double r_ij = 0;

                    for (int coord = 0; coord < 3; coord++){
                        
                        d_xyz[coord] = 2*(pos[3 * n + coord] - pos[3 * i_neig + coord]); // vector ij' twice the lenght of classical ij vector
                        r_ij += d_xyz[coord]*d_xyz[coord];
                    }

                    r_ij = sqrt(r_ij);

                    if (r_ij < geomParams.kappa*geomParams.h){ // if new imaginary j' particle in support domain of i

                        vector<double> new_gradW(3);
                        double deriv = derive_wendland_quintic(r_ij, geomParams, simParams);
                    
                        for (int coord = 0; coord < 3; coord++){
                            new_gradW[coord] = (d_xyz[coord] / r_ij) * deriv;
                            normal[3*n+coord] -= (0 - color[n])*(mass[n]/rho[n])*(new_gradW[coord]);
                            //normal[3*n+coord] -= (0 - color[n])*(mass[i_neig]/rho[i_neig])*(1);

                        }
                    }
                }
            }     
            
                    
        }

        
        double norm = 0;
       // cout << "n = " << n << endl;
        for (int coord = 0; coord < 3; coord++){
            norm += normal[3*n+coord]*normal[3*n+coord];
        }

        norm = sqrt(norm);
        

        R[n] = (norm <= 0.01/geomParams.h)? 0: 1;
        cout << "norm = " << norm << endl;

        // normalize normal 

        //cout << "n : "<< n << endl;
        //cout << "norm= " << norm << endl;

        for (int coord = 0; coord < 3; coord++){
                  
            //cout << "before normalize" << endl;
            //cout << "normal[3*n + " << coord << "]= " << normal[3*n+coord] << endl;


            normal[3*n+coord] /= norm;
            //cout << "after normalize" << endl;
            //cout << "normal[3*n + " << coord << "]= " << normal[3*n+coord] << endl;
            //cout << "\n";

        }
        
    }

    //printMatrix(gradW_matrix, gradW_matrix.size(), "gradW_matrix");
    //printArray(normal, normal.size(), "normal"); 
    //printArray(R, R.size(), "R"); 

  
    double Kappa = 0, L = 0;

    // Calculation of curvature
    //#pragma omp parallel for
    for (int n = 0; n < simParams.nb_moving_part; n++){
    //for (int n = 0; n < 40; n++){

        int size_neighbours = nb_neighbours[n];
        //cout << "---------------------"<< endl;
        //cout << "n : "<< n << endl;

        for(int idx = 0; idx < size_neighbours; idx++){

            int i_neig = neighbours[100*n + idx];
          //  cout << "i_neig : "<< i_neig << endl;

            // classic curvature calculation  
            double dot_product = 0;
            for (int coord = 0; coord < 3; coord++){

             //   cout << "coord = " << coord << endl;
             //cout << "normal[3*i_neig+coord] - normal[3*n+coord] = " << normal[3*i_neig+coord] - normal[3*n+coord] << endl;
            // cout << "gradW[n][3*idx+coord] = " << gradW[n][3*idx+coord] << endl;
                
                dot_product += (normal[3*i_neig + coord] - normal[3*n + coord])*gradW[n][3*idx + coord];
            }

            //cout << "i_neig : "<< i_neig  << " (" << idx << "/"  << size_neighbours << ")" << endl;

            //cout << "dot_product = " << dot_product << endl;
            //cout << "R[n] :" << R[n] << endl;
            //cout << "R[i_neig] :" << R[i_neig] << endl;
            //cout << "mass[i_neig] :" << mass[i_neig] << endl;
            //cout << "rho[i_neig] :" << rho[i_neig] << endl;

            double a = R[n]*R[i_neig];
            double b = mass[i_neig]/rho[i_neig];
            double c = W[n][idx];

            //cout << "a : " << R[n]*R[i_neig] << endl;
            //cout << "b : " << mass[i_neig]/rho[i_neig] << endl;
            //cout << "c : " << dot_product << endl;
            //cout << "c : " << W[n][idx] << endl;
            //cout << "a*b*c : " << a*b*c << endl;


            Kappa += (R[n]*R[i_neig])*(mass[i_neig]/rho[i_neig])*dot_product;
            //cout << "Kappa  = " << Kappa << endl;
            //cout << "\n";

            L += (R[n]*R[i_neig])*(mass[i_neig]/rho[i_neig])*W[n][idx];
            //cout << "L = " << L << endl;

            
            // imaginary particle contribution
            if (track_particle[n] == 1 && track_particle[i_neig] == 0){ // particle on free surface but not neighbour 

                vector<double> new_normal(3);
                double norm = 0;

                for (int coord = 0; coord < 3; coord++){
                    new_normal[coord] = 2*normal[3*n+coord] - normal[3*i_neig+coord]; 
                    norm += new_normal[coord]*new_normal[coord];
                }

                // normalize new vector
                for (int coord = 0; coord < 3; coord++)
                    new_normal[coord] = new_normal[coord]/norm;

                // contribution of normal ij'
                double dot_product = 0;

                for (int coord = 0; coord < 3; coord++)
                    dot_product += (new_normal[coord] - normal[3*n+coord])*gradW[n][3*idx+coord];

                Kappa += (R[n]*R[i_neig])*(mass[i_neig]/rho[i_neig])*dot_product;
                cout << "Second operation on kappa :" << Kappa << endl;


                

                L += (R[n]*R[i_neig])*(mass[i_neig]/rho[i_neig])*W[n][idx];
                
                continue;   
            }

            if (N[n]){ // at least one free surface neighbour in support domain
                if (track_particle[n] == 0 && track_particle[i_neig] == 1){ // neighbour on free surface but not particle 

                    vector<double> new_normal(3);
                    double norm = 0;

                    for (int coord = 0; coord < 3; coord++){
                        new_normal[coord] = 2*normal[3*i_neig+coord] - normal[3*n+coord]; 
                        norm += new_normal[coord]*new_normal[coord];
                    }

                    // normalize new vector
                    for (int coord = 0; coord < 3; coord++)
                        new_normal[coord] = new_normal[coord]/norm;

                    // contribution of normal ij'
                    double dot_product = 0;

                    for (int coord = 0; coord < 3; coord++)
                        dot_product += (new_normal[coord] - normal[3*n+coord])*gradW[n][3*idx+coord];

                    Kappa += (R[n]*R[i_neig])*(mass[i_neig]/rho[i_neig])*dot_product;
                    cout << "Third operation on kappa :" << Kappa << endl;


                    L += (R[n]*R[i_neig])*(mass[i_neig]/rho[i_neig])*W[n][idx];

                }
            }              
        }

       

        //cout << "n = " << n << endl;
        //cout << "Kappa (fin) = " << Kappa << endl;
        //cout << "L = " << L << endl;
        //cout << "\n";

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



