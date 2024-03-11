#include <stdio.h>
#include <string.h>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <cassert>

#include "generate_particle.h"
#include "sorted_list.h"
#include "gradient.h"
#include "initialize.h"
#include "export.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "nlohmann/json.hpp"
using json = nlohmann::json;
using namespace std;



/**
 * @brief Dummy SPH simulation testing paraview export formats
 */

int main(int argc, char *argv[]){

    /*---------------------------- SETTING COMPILATION PARAMETERS -----------------------------------*/

    #ifdef _OPENMP
        std::cout << "OpenMP available: OMP_NUM_THREADS=" << omp_get_max_threads() << "\n";
    #else
        std::cout << "OpenMP not available.\n";
    #endif

    #ifdef NDEBUG
        // code has been configured with "cmake -DCMAKE_BUILD_TYPE=Release .."
        std::cout << "code built in RELEASE mode.\n";
    #else
        // code has been configured with "cmake .."
        std::cout << "code built in DEBUG mode.\n";
    #endif

    if (argc != 2)
    {
        std::cout << "\nusage: " << argv[0] << " <param.json>\n\n";
        return EXIT_FAILURE;
    }

    /*---------------------------- INPUT PARAMETERS FROM JSON FILES -----------------------------------*/

    std::ifstream inputf(argv[1]);
    json data = json::parse(inputf);

    std::cout << argv[1] << ":\n" << data.dump(4) << std::endl;     // Print input data to screen

    std::vector<double> o = data["o"];
    std::vector<double> L = data["L"];
    double s = data["s"];
    int nstepT = data["nstepT"];
    std::vector<double> o_d = data["o_d"];
    std::vector<double> L_d = data["L_d"];
    
    std::string state_equation_chosen;

    for (auto& it : data["stateEquation"].items()){
        if (it.value    () == 1){
            state_equation_chosen = it.key();
        }
    }
    
    cout << "state equation chose : " << state_equation_chosen << " \n" <<endl; 

    unsigned Nx, Ny, Nz ;    // Number of cells (in each direction)

    int kappa = data["kappa"];
    double h = 1.2*s;
    double R = 8.314; // [J/(K.mol)]
    double g = -9.81; // [m/s²]

    Nx = (int) L_d[0]/(kappa*h);
    Ny = (int) L_d[1]/(kappa*h);
    Nz = (int) L_d[2]/(kappa*h);
    cout << " kappa * h =" << kappa*h << endl;
    printf("(Nx, Ny, Nz) = (%d, %d, %d) \n", Nx,Ny,Nz);

    vector<double> pos_arr;  
    vector<vector<unsigned>> cell_pos(Nx*Ny*Nz);
    meshcube(&o[0], &L[0], s, pos_arr); // Initialise particles in the domain

    unsigned nb_particles = pos_arr.size()/3;

    double rho_init = data["rho"];
    vector<double> u_init = data["u"];
    vector<double> mass_arr(nb_particles), u_arr(3*nb_particles), drhodt_arr(nb_particles), rho_arr(nb_particles), dudt_arr(3*nb_particles), p_arr(nb_particles);
    vector<vector<double>> gradW_matrix, artificial_visc_matrix(nb_particles); 
    vector<vector<unsigned>>  neighbours_matrix(nb_particles);

    double dt = 0.1;
    initializeRho(rho_arr,rho_init);
    initializeMass(rho_arr, s, mass_arr);
    initializeVelocity(u_arr, u_init);

    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;
    
    vectors["position"] = &pos_arr;
    vectors["velocity"] = &u_arr;

    double rho_0 = 1000, c_0 = 340.0, T = 273, M = 200, gamma = 7;

for(size_t i = 0; i<artificial_visc_matrix.size();i++){
    for(size_t j=0 ; artificial_visc_matrix[i].size();j++){
        artificial_visc_matrix[i][j] = 0.0;
    }


}  

    
  /*---------------------------- SPH ALGORITHM  -----------------------------------*/

    vector<double> u_temp(3*nb_particles), rho_temp(nb_particles), pos_temp(nb_particles);

    for (int t = 0; t < nstepT; t++){
        
        // Apply the linked-list algorithm
        findNeighbours(pos_arr, cell_pos, neighbours_matrix, &L_d[0], Nx, Ny, Nz, h, kappa);

        // Compute ∇_a(W_ab) for all particles
        gradW(pos_arr, neighbours_matrix, gradW_matrix, h, Nx, Ny, Nz);

        // Compute D(rho)/Dt for all particles
        continuityEquation(pos_arr ,u_arr, neighbours_matrix, gradW_matrix, drhodt_arr, rho_arr, mass_arr, h);

        // Compute D(u)/Dt for all particles
        momentumEquation(neighbours_matrix, mass_arr, gradW_matrix, dudt_arr, artificial_visc_matrix, rho_arr, rho_0, c_0, p_arr, R, T, M, gamma, state_equation_chosen);

        // Update density, velocity and position for each particle (Euler explicit scheme)
        for(size_t pos = 0; pos < nb_particles; pos++ ){

            rho_arr[pos] = rho_arr[pos] + dt*drhodt_arr[pos];
            
            //cout << " drho dt est de " << drhodt_arr[pos]<< endl;

            for (size_t cord = 0; cord < 2; cord++){
                pos_arr[3*pos+cord] = pos_arr[3*pos+cord] + dt*u_arr[3*pos+cord];
                
                u_arr[3*pos+cord] = u_arr[3*pos+cord] + dt*dudt_arr[3*pos+cord];
                
                u_arr[3*pos+2] = u_arr[3*pos+2] + dt*g*1000;
                pos_arr[3*pos+cord] = (pos_arr[3*pos+cord] < 0.0) ? 0.0 : pos_arr[3*pos+cord];
                pos_arr[3*pos+cord] = (pos_arr[3*pos+cord] > L_d[cord]) ? L_d[cord] : pos_arr[3*pos+cord];
                u_arr[3*pos+cord] = (pos_arr[3*pos+cord] < 0.0) ? 0.0 : u_arr[3*pos+cord];
                u_arr[3*pos+cord] = (pos_arr[3*pos+cord] > L_d[cord]) ? 0.0 : u_arr[3*pos+cord];
            }
            

        }
         
        for(size_t i = 0 ; i<neighbours_matrix.size(); i++ ){
            
            neighbours_matrix[i].clear();
            gradW_matrix[i].clear();
        }
        gradW_matrix.clear();
        for(size_t i = 0 ; i<cell_pos.size(); i++ ){
            
            cell_pos[i].clear();
        }
        export_particles("sph", t, pos_arr, scalars, vectors);
    }
}


