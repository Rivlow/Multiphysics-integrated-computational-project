#include <stdio.h>
#include <string.h>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>

#include "generate_particle.h"
#include "sorted_list.h"
#include "gradient.h"
#include "initialize.h"
#include "export.h"
#include "Kernel_functions.h"
#include "tools.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#include <chrono>

#include "nlohmann/json.hpp"
using json = nlohmann::json;
using namespace std;


int main(int argc, char *argv[])
{

    /*------------- SETTING COMPILATION PARAMETERS/Folders needed ---------------------------*/
auto t0 = std::chrono::high_resolution_clock::now();

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


    /*---------------------------- INPUT PARAMETERS FROM JSON FILES ----------*/

    std::ifstream inputf(argv[1]);
    json data = json::parse(inputf);

    std::cout << argv[1] << ":\n"
              << data.dump(4) << std::endl; // Print input data to screen


    std::string state_equation_chosen; /// p expl "Ideal gaz law"
    std::string state_initial_condition;

    for (auto &it : data["stateEquation"].items())
    {
        if (it.value() == true)
        {
            state_equation_chosen = it.key();
        }
    }

    for (auto &it : data["initialCondition"].items())
    {
        if (it.value() == true)
        {
            state_initial_condition = it.key();
        }
    }

    createOutputFolder();
    clearOutputFiles();


    SimulationData params;
    params.kappa = data["kappa"];
    params.dt = data["dt"];
    params.s = data["s"];
    params.h = 1.2 * params.s; 
    params.nstepT = data["nstepT"]; 
    params.o = data["o"].get<vector<double>>();
    params.L = data["L"].get<vector<double>>();
    params.o_d = data["o_d"].get<vector<double>>();
    params.L_d = data["L_d"].get<vector<double>>();
    params.u_init = data["u_init"].get<vector<double>>();
    params.alpha = data["alpha"];
    params.beta = data["beta"];
    params.c_0 = data["c_0"];
    params.rho_moving = data["rho_moving"];
    params.rho_fixed = data["rho_fixed"];
    params.rho_0 = data["rho_0"];
    params.M = data["M"];
    params.T = data["T"];
    params.gamma = data["gamma"];
    params.nsave = data["nsave"];
    params.state_equation = state_equation_chosen; 
    params.PRINT = data["print_debug"];
    std::string state_initial_condition;
    
  

    /*---------------------------- INITIALIZATION OF VARIABLES USED ----------*/



    // Constants
                 
   

    // Number of cells 
    int Nx = params.L_d[0] / (params.kappa * params.h);
    int Ny = params.L_d[1] / (params.kappa * params.h);
    int Nz = params.L_d[2] / (params.kappa * params.h);

    cout << "state equation chosen : " << state_equation_chosen << " \n" << endl;
    cout << " kappa * h =" << params.kappa * params.h << endl;
    printf("(Nx, Ny, Nz) = (%d, %d, %d) \n", Nx, Ny, Nz);

    // Nb of particles along each direction from target size "s"
    size_t nb_particles = evaluateNumberParticles(params);

    // Vector used (and labelled w.r.t particles location)
    vector<double> pos_array, type_arr;
    vector<vector<int>> cell_matrix(Nx * Ny * Nz);

    // Variables defined to used "export.cpp"
    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;

    /*---------------------------- SPH ALGORITHM  ----------------------------*/

    // Initialization of the problem (moving particles and fixed particles)
    meshcube(params, pos_array, type_arr);
    size_t nb_moving_part = pos_array.size() / 3;
    meshBoundary(params, pos_array, type_arr);

    int nb_tot_part = pos_array.size()/3;

    std::cout << "nb_moving_part = " << nb_moving_part << std::endl;
    std::cout << "nb_tot_part = " << nb_tot_part << std::endl;
    std::cout << "s=" << params.s << std::endl;
    std::cout << "kappa=" << params.kappa << std::endl;
    std::cout << "h=" << params.h << std::endl;

    vector<double> mass_array(nb_tot_part),
                   u_array(3 * nb_tot_part),
                   drhodt_array(nb_tot_part),
                   rho_array(nb_tot_part),
                   dudt_array(3 * nb_tot_part, 0.0),
                   p_array(nb_tot_part),
                   c_array(nb_tot_part),
                   grad_sum(nb_tot_part);

    vector<vector<double>> artificial_visc_matrix(nb_tot_part),
                           gradW_matrix(nb_tot_part);
    vector<vector<int>> neighbours_matrix(nb_tot_part);
    std::vector<double> nvoisins(nb_tot_part, 0.0); //RB


    scalars["type"] = &type_arr;
    scalars["mass_array"] = &mass_array;
    scalars["rho_array"] = &rho_array;
    scalars["p_array"] = &p_array;
    scalars["drhodt_array"] = &drhodt_array;
    scalars["nvoisins"] = &nvoisins;
    scalars["grad_sum"] = &grad_sum;
    vectors["position"] = &pos_array;
    vectors["u_array"] = &u_array;
    vectors["dudt_array"] = &dudt_array;


    initializeRho(params, pos_array, rho_array,
                  nb_moving_part,
                  state_equation_chosen, state_initial_condition);

    initializeMass(params, rho_array, mass_array);
    initializeVelocity(params, u_array, nb_moving_part);
    initializeViscosity(params, artificial_visc_matrix);


    for (int t = 0; t < nstepT; t++)
    {

        // compute dt_max = min_a h_a/F
        // double dt_max = sqrt(kappa * h / abs(g)); // CFL condition
        // std::cout << "dt_max = " << dt_max << "    dt = " << dt << std::endl;

        // Apply the linked-list algorithm
        findNeighbours(cell_matrix, neighbours_matrix,
                       pos_array, L_d,
                       nb_moving_part, Nx, Ny, Nz, h, kappa, 
                       PRINT); 


        // Compute ∇_a(W_ab) for all particles
        gradW(gradW_matrix, neighbours_matrix, 
              pos_array,
              nb_moving_part, h, Nx, Ny, Nz,
              PRINT); 


        // Compute D(rho)/Dt for all particles
        continuityEquation(neighbours_matrix, gradW_matrix, 
                           pos_array, u_array, drhodt_array, rho_array, mass_array,
                           nb_moving_part, h,
                           PRINT); 


        // Compute D(u)/Dt for all particles
        momentumEquation(neighbours_matrix, gradW_matrix, artificial_visc_matrix,
                         mass_array, dudt_array, rho_array, p_array, c_array, pos_array, u_array,
                         nb_moving_part,
                         rho_0, c_0, gamma, R, T, M, g, t, alpha, beta, h,
                         state_equation_chosen, PRINT); 
        

        // Update density, velocity and position for each particle (Euler explicit scheme)
        update(pos_array, u_array, rho_array, drhodt_array, dudt_array, dt, PRINT);
        
        // Clear matrices and reset to 0 arrays
        clearAllVectors(artificial_visc_matrix, neighbours_matrix, cell_matrix, gradW_matrix,
                        drhodt_array, dudt_array,
                        PRINT);

        if(t % nsave == 0){
            export_particles("../../output/sph", t, pos_array, scalars, vectors);
        }

    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    std::cout << "duration: " << delta_t << "s.\n";

    std::cout << "\n Simulation done." << std::endl;
    return 0;
}
