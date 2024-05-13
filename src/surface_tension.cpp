#include <stdio.h>
#include <vector>
#include <cmath>
#include <omp.h>

#include "gradient.h"
#include "find_neighbours.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"
#include "surface_tension.h"
#include "Kernel.h"

using namespace std;



void surfaceTension(SimulationData& simParams,
                    GeomData &geomParams,
                    ThermoData &thermoParam,
                    vector<double> nb_neighbours,
                    vector<int> neighbours,
                    vector<int> &track_surface,
                    vector<double> &N_smoothed,
                    vector<vector<double>> gradW_matrix,
                    vector<vector<double>> W_matrix,
                    vector<double> mass,
                    vector<double> rho,
                    vector<double> pos,
                    vector<double> &F_vol,
                    vector<double> type){


    vector<double> normal(3*simParams.nb_moving_part,0.0);
    vector<double> div_r(3*simParams.nb_moving_part,0.0);
    vector<double> c(simParams.nb_moving_part,0.0);
    vector<double> R(simParams.nb_moving_part,0.0);
    vector<double> N(simParams.nb_moving_part,0.0);


    cout << "first loop " << endl;
    // First loop compute the scalar div(r) 
    // and determine the "N" for each particle (N = 1 indicates particle at free surface, otherwise N = 0)
    // The free surface is detected by both math (first) and geom method (second)
    //printArray(track_surface, track_surface.size(), "tracksur");
    int nb_one = 0;
    int nb_zero = 0;
    
    for(int n = 0; n < simParams.nb_moving_part; n++){

        vector<double> &gradW = gradW_matrix[n];
        int size_neighbours = nb_neighbours[n];
        
        for(int idx = 0; idx < size_neighbours; idx++){ 
            
            int i_neig = neighbours[100*n + idx]; 
            double m_j = mass[i_neig];
            double rho_j = rho[i_neig];
            double dot_product = 0; //r_ij*gradW_ij
            
            for(int coord = 0; coord < 3; coord ++){
                
                double r_ij = pos[3*n+coord] - pos[3*i_neig+coord];
                
                double gradW_ij = gradW[3*idx+coord];
                
                dot_product += r_ij*m_j*gradW_ij/rho_j;
                //normal[3*n+coord] += r_ij*m_j*gradW_ij/rho_j;
                
            }
            //cout << rho_j << endl;
           // cout<< " dot_prod : " << dot_product << endl;
            div_r[n] = dot_product;
            //cout<< " apres dot prot "  << endl;
        }
        //cout<<"après 2eme " << endl;
        // mathematical check
        if (1){ // particle n is potentially located on free surface

            int count = 0; // geom check

            for (int i = 0; i < 8; i++){

                if (track_surface[8*n + i] == 0)
                    count++;   
            }

            if (count > 0)
            {
                N[n] = 1;
                nb_one++;
            }

        }
        nb_zero++;
        
    }
    //printArray(track_surface,track_surface.size(), "track_surface");
    
    cout << "nb_one " << nb_one  << " nb_zero " << nb_zero<< endl;

    // Second loop evaluated the smoothed color function "c"
    #pragma omp parallel for 
    for(int n = 0; n < simParams.nb_moving_part; n++){

        int size_neighbours = nb_neighbours[n];
        double sum = 0;
        
        if (N[n] == 1){
            for(int idx = 0; idx < size_neighbours; idx++){ 
                
                int i_neig = neighbours[100*n + idx];
                double W_ij = W_matrix[n][idx];
                double m_j = mass[i_neig];
                double rho_j = rho[i_neig];

                sum += m_j*W_ij/rho_j;
            }
            
            c[n] = sum;
        }

        else
            c[n] = 1;
    }
   
    //cout << "Second loop ok " << endl;


    // Third loop evaluated the normal vector "n"
    #pragma omp parallel for 
    for(int n = 0; n < simParams.nb_moving_part; n++){

        int size_neighbours = nb_neighbours[n];

        for(int idx = 0; idx < size_neighbours; idx++){ 

            int i_neig = neighbours[100*n + idx];
           
            double m_j = mass[i_neig];
            double rho_j = rho[i_neig];

            for (int coord = 0; coord < 3; coord++)
                normal[3*n + coord] -= m_j*gradW_matrix[n][3*idx+coord]*(c[i_neig] - c[n])/rho_j;
        }
    }
    
    //cout << "Third loop ok " << endl;
    //printArray(N,N.size(),"N");
    // Fourth loop evaluate the surface curvature
    //#pragma omp parallel for 
    for(int n = 0; n < simParams.nb_moving_part; n++){

        int size_neighbours = nb_neighbours[n];

        // Evaluate normal n
        for(int idx = 0; idx < size_neighbours; idx++){ 

            int i_neig = neighbours[100*n + idx];
           
            double m_j = mass[i_neig];
            double rho_j = rho[i_neig];

            for (int coord = 0; coord < 3; coord++){
                normal[3*n + coord] -= m_j*gradW_matrix[n][3*idx+coord]*(c[i_neig] - c[n])/rho_j;
                //cout << " normal[3*n + coord] : "<< normal[3*n + coord] << "    ";
            }
                //cout << endl;
        }

        double norm_i = 0, norm_j = 0;
        double k_i = 0, L_i = 0;
        
        // Evaluate surface curvature
        for(int idx = 0; idx < size_neighbours; idx++){ 

            int i_neig = neighbours[100*n + idx];
            double m_j = mass[i_neig];
            double rho_j = rho[i_neig];
            double N_i = N[n];
            double N_j = N[i_neig];

            for (int coord = 0; coord < 3; coord++){
                norm_i +=normal[3*n + coord]*normal[3*n + coord];
                norm_j +=normal[3*i_neig + coord]*normal[3*i_neig + coord];
            }


            norm_i= sqrt(norm_i);
            norm_j= sqrt(norm_j);
            //cout << "norm_i : "<< norm_i << ", norm_j : " << norm_j << endl;
            // find particles with reliable normal vectors (truncate solution using R_i)
            double R_i = (norm_i > 0.01/geomParams.h)? 1 : 0;

            // If particle i or j at the surface, create imaginary particle 
            vector<double> n_imag(3);
            double norm_imag = 0;
            if(norm_i != 0 && norm_j !=0){

                if (N_i == 1.0 && N_j == 0.0){ // (particle i at surface and j in the fluid)
                //cout << "test1" << endl;
                    for (int coord = 0; coord < 3; coord++){
                        n_imag[coord] = 2*normal[3*n + coord]/norm_i - normal[3*i_neig + coord]/norm_j;
                        norm_imag += n_imag[coord]*n_imag[coord];
                    }

                    norm_imag = sqrt(norm_imag);
                    double R_imag = (norm_imag > 0.01/geomParams.h)? 1 : 0;
                    double dot_product = 0;

                    for (int coord = 0; coord < 3; coord++){
                        double delta_norm = n_imag[coord]/norm_imag - normal[3*n+coord]/norm_i;
                        dot_product += delta_norm*(-1*gradW_matrix[n][3*idx+coord]);
                    }

                    k_i -= min(R_i, R_imag)*dot_product*m_j/rho_j;
                    L_i += min(R_i, R_imag)*m_j*W_matrix[n][idx]/rho_j; // kernel gradient correction to counteract truncated solution
                }

                else if (N_j == 1.0 && N_i == 0.0){ 
                    continue;// (particle j at surface and i in the fluid)
                    //cout << "ahhhhh1" << endl;
                    for (int coord = 0; coord < 3; coord++){
                        n_imag[coord] = 2*normal[3*i_neig + coord]/norm_j - normal[3*n + coord]/norm_i;
                        norm_imag += n_imag[coord]*n_imag[coord];
                    }

                    norm_imag = sqrt(norm_imag);
                    double R_imag = (norm_imag > 0.01/geomParams.h)? 1 : 0;
                    double dot_product = 0;

                    for (int coord = 0; coord < 3; coord++){
                        double delta_norm = n_imag[coord]/norm_imag - normal[3*n+coord]/norm_i;
                        dot_product += delta_norm*(-1*gradW_matrix[n][3*idx+coord]);// a chané ça car grad pas bon
                    }

                    k_i -= min(R_i, R_imag)*dot_product*m_j/rho_j;
                    L_i += min(R_i, R_imag)*m_j*W_matrix[n][idx]/rho_j; // kernel gradient correction to counteract truncated solution

                }
                else{ // neither particle i,j on the surface
                    
                    double R_j = (norm_j > 0.01/geomParams.h)? 1 : 0;           
                    double dot_product = 0;
                    
                    for (int coord = 0; coord < 3; coord++){
                        double delta_norm = normal[3*i_neig + coord]/norm_j - normal[3*n+coord]/norm_i;
                        dot_product += delta_norm*(gradW_matrix[n][3*idx+coord]);
                        /*if(abs(dot_product)<0.000001 && abs(dot_product)>0){
                            cout << " R_j  " << R_j <<"  norm_j  " << norm_j <<  endl;
                        cout<<"dot product : " << dot_product << " delta_norm : " << delta_norm << " gradW_matrix[n][3*idx+coord] " << gradW_matrix[n][3*idx+coord]<<endl;
                        cout<<"n " << n << " ineig " << i_neig << endl;
                        cout << "normal[3*i_neig + coord]/norm_j : "  << normal[3*i_neig + coord] << " normal[3*n+coord]/norm_i " << normal[3*n+coord]/norm_i << endl;
                        }*/

                    }
                    k_i -= R_i*R_j*dot_product*m_j/rho_j;
                    L_i += R_i*R_j*m_j*W_matrix[n][idx]/rho_j; 
                    
                    cout << "   n   "<< n<<  "      dot_product   " << dot_product << "   k_i     "<< k_i << "    R_i     " << R_i << "      R_j     " << R_j << "    m_j    " <<  m_j << endl;// kernel gradient correction to counteract truncated solution
                }
            }
            if(L_i){
                k_i  = k_i/L_i;
                
            }   
            
            
            
            

            

            
        }
        double F_res = 0;
        for (int coord = 0; coord < 3; coord++){
                    F_vol[3*n + coord] += thermoParam.sigma*k_i*normal[3*n + coord]/rho[n];
                    //cout<<" k_i " << k_i << "   normal[3*n + coord]     " << normal[3*n + coord] <<  "    rho[n]    "  << rho[n] << endl;
                    F_res += F_vol[3*n + coord]*F_vol[3*n + coord];
                }
        simParams.F_st_max = sqrt(F_res);
            
    }  
    //cout << "Fourth loop ok " << endl;

    
}

/*double W_coh(double r, double h){
    double W = 0.0;
    double cst = 32/(M_PI*h*h*h*h*h*h*h*h*h);
    if(2.0*r>h && r<=h){
        W = cst*(h-r)*(h-r)*(h-r)*r*r*r;
    }
    else if(r>0 && 2.0*r<=h){
        W = cst *( 2.0*(h-r)*(h-r)*(h-r)*r*r*r - h*h*h*h*h*h/64);
    }
    else{
        W = 0.0;
    }
    return W;
}

void surfaceTension(SimulationData& simParams,
                    GeomData &geomParams,
                    ThermoData &thermoParam,
                    vector<double> nb_neighbours,
                    vector<int> neighbours,
                    vector<vector<double>> gradW_matrix,
                    vector<double> mass,
                    vector<double> rho,
                    vector<double> pos,
                    vector<double> &F_vol){

    vector<double> normal(3*simParams.nb_moving_part,0.0);
    #pragma omp parallel for 
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
    #pragma omp parallel for 
    for(int n = 0; n<simParams.nb_moving_part; n++){
    
        int size_neighbours = nb_neighbours[n];
        
        for(int idx = 0; idx<size_neighbours; idx++){
            if(type[neighbours[100*n + idx]] == 1){
                int i_neig = neighbours[100*n + idx];
                double K_ij = 2*thermoParam.rho_0/(rho[n]+rho[i_neig]);
                double r_ab = 0;
                vector<double> d_xyz(3);
                
                for (int coord = 0; coord < 3; coord++){
                    
                    d_xyz[coord] = pos[3 * n + coord] - pos[3 * i_neig + coord];
                    r_ab += d_xyz[coord]*d_xyz[coord];
                }
                
                r_ab = sqrt(r_ab);
                double W_ab = W_coh(r_ab,geomParams.kappa*geomParams.h);
                double m_a = mass[n];
                double m_b = mass[i_neig];
                double F_res = 0;
                
                for (int coord = 0; coord < 3; coord++){
                    
                    F_vol[3*n + coord] += -K_ij*(alpha * m_a * m_b * d_xyz[coord]*W_ab/r_ab 
                                    + alpha*(normal[3*n+coord]-normal[3*i_neig+coord]));
                
                    F_res += F_vol[3*n + coord]*F_vol[3*n + coord];
                }
                
                simParams.F_st_max = sqrt(F_res);
            }
        }
            
    }

    












    
*/