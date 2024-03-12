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
#include "Kernel_functions.h"


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

    std::cout << argv[1] << ":\n" << data.dump(4) << std::endl; // Print input data to screen

    std::vector<double> o = data["o"];
    std::vector<double> L = data["L"];
    double s = data["s"];
    int nstepT = data["nstepT"];
    std::vector<double> o_d = data["o_d"];
    std::vector<double> L_d = data["L_d"];
    
    std::string state_equation_chosen, state_initial_condition;

    for (auto& it : data["stateEquation"].items()){
        if (it.value    () == 1){
            state_equation_chosen = it.key();
        }
    }

    for (auto& it : data["initialCondition"].items()){
        if (it.value    () == 1){
            state_initial_condition = it.key();
        }
    }
    
    cout << "state equation chose : " << state_equation_chosen << " \n" <<endl; 

    int kappa = data["kappa"];
    double alpha = data["alpha"];
    double beta = data["beta"];
    double c_0 = data["c_0"];

    double rho_init = data["rho"];
    double rho_0 = data["rho_0"];
    double M = data["M"];
    double T = data["T"];
    double gamma = data["gamma"];
    vector<double> u_init = data["u"];

    /*---------------------------- INITIALIZATION OF VARIABLES USED -----------------------------------*/

    // Constants
    double h = 1.2*s;
    double R = 8.314; // [J/(K.mol)]
    double g = 9.81; // [m/s²]
    double dt = 0.005;

    // Number of cells (in each direction) 
    unsigned Nx, Ny, Nz ;    
    Nx = (int) L_d[0]/(kappa*h);
    Ny = (int) L_d[1]/(kappa*h);
    Nz = (int) L_d[2]/(kappa*h);

    cout << " kappa * h =" << kappa*h << endl;
    printf("(Nx, Ny, Nz) = (%d, %d, %d) \n", Nx,Ny,Nz);

    // Nb of particles along each direction from target size "s"
    int nb_particles = evaluateNumberParticles(L, s);

    // Vector used (and labelled w.r.t particles location)
    vector<double> pos_arr;  
    vector<vector<unsigned>> cell_matrix(Nx*Ny*Nz);
    vector<double> mass_arr(nb_particles), u_arr(3*nb_particles), drhodt_arr(nb_particles), rho_arr(nb_particles), dudt_arr(3*nb_particles), p_arr(nb_particles);
    vector<vector<double>> gradW_matrix, artificial_visc_matrix(nb_particles); 
    vector<vector<unsigned>>  neighbours_matrix(nb_particles), neighbours_matrix_1(nb_particles);

    // Variables defined to used "export.cpp"
    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;
    vectors["position"] = &pos_arr;
    vectors["velocity"] = &u_arr;

  /*---------------------------- SPH ALGORITHM  -----------------------------------*/

    // Initialization of the problem
    meshcube(o, L, s, pos_arr); 
    initializeRho(pos_arr, rho_arr, rho_init, rho_0, c_0, M, g, R, T, gamma, state_equation_chosen, state_initial_condition);
    initializeMass(rho_arr, s, mass_arr);
    initializeVelocity(u_arr, u_init);
    initializeViscosity(artificial_visc_matrix);

    for (int t = 0; t < nstepT; t++){

        findNeighbours(pos_arr, cell_matrix, neighbours_matrix, L_d, Nx, Ny, Nz, h, kappa); // Apply the linked-list algorithm
        gradW(pos_arr, neighbours_matrix, gradW_matrix, h, Nx, Ny, Nz); // Compute ∇_a(W_ab) for all particles
        continuityEquation(pos_arr ,u_arr, neighbours_matrix, gradW_matrix, drhodt_arr, rho_arr, mass_arr, h); // Compute D(rho)/Dt for all particles
        setPressure(p_arr, rho_arr, rho_0, c_0, R, T, M, gamma, state_equation_chosen); // Compute pressure for all particles
        setArtificialViscosity(t, artificial_visc_matrix, pos_arr, neighbours_matrix, rho_arr, u_arr,
                               alpha, beta, rho_0, c_0, gamma, R, T, M, h, state_equation_chosen); // Compute artificial viscosity Π_ab for all particles
        momentumEquation(neighbours_matrix, mass_arr, gradW_matrix, dudt_arr, artificial_visc_matrix, rho_arr, 
                        rho_0, c_0, p_arr, R, T, M, gamma, g, state_equation_chosen); // Compute D(u)/Dt for all particles

        // Update density, velocity and position for each particle (Euler explicit scheme)
        for(size_t pos = 0; pos < nb_particles; pos++ ){

            rho_arr[pos] = rho_arr[pos] + dt*drhodt_arr[pos];

            for (size_t cord = 0; cord < 3; cord++){

                pos_arr[3*pos+cord] += dt*u_arr[3*pos+cord];
                u_arr[3*pos+cord] += dt*dudt_arr[3*pos+cord];

                // Check boundaries (temporary)
                u_arr[3*pos+cord] = (pos_arr[3*pos+cord] < 0.0) ? 0.0 : u_arr[3*pos+cord];
                u_arr[3*pos+cord] = (pos_arr[3*pos+cord] > L_d[cord]) ? 0.0 : u_arr[3*pos+cord];
                pos_arr[3*pos+cord] = (pos_arr[3*pos+cord] < 0.0) ? 0.0 : pos_arr[3*pos+cord];
                pos_arr[3*pos+cord] = (pos_arr[3*pos+cord] > L_d[cord]) ? L_d[cord] : pos_arr[3*pos+cord];
            }
        }

        clearAllVectors(artificial_visc_matrix, neighbours_matrix, cell_matrix, gradW_matrix);
        export_particles("sph", t, pos_arr, scalars, vectors);
    }
}


