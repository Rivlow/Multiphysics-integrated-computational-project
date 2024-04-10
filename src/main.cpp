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

    std::vector<double> o = data["o"];
    std::vector<double> L = data["L"];
    double s = data["s"];
    int nstepT = data["nstepT"];
    std::vector<double> o_d = data["o_d"];
    std::vector<double> L_d = data["L_d"];

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

    cout << "state equation chosen : " << state_equation_chosen << " \n"
         << endl;

    int kappa = data["kappa"];
    double alpha = data["alpha"];
    double beta = data["beta"];
    double c_0 = data["c_0"];
    double rho_init = data["rho_moving"];
    double rho_fixed = data["rho_fixed"];
    double rho_0 = data["rho_0"];
    double M = data["M"];
    double T = data["T"];
    double gamma = data["gamma"];
    double dt = data["dt"]; // RB
    int nsave = data["nsave"]; // RB
    vector<double> u_init = data["u"];         // RB: inutile

    /*---------------------------- INITIALIZATION OF VARIABLES USED ----------*/

    // Debug variable
    const bool PRINT = data["print_debug"];

    // Constants
    double h = 1.2 * s;                    
    double R = 8.314; // [J/(K.mol)]
    double g = -9.81; // [m/s²]

     SimulationData params = {

    cout << " kappa * h =" << kappa * h << endl;
    printf("(Nx, Ny, Nz) = (%d, %d, %d) \n", Nx, Ny, Nz);

    // Nb of particles along each direction from target size "s"
    size_t nb_particles = evaluateNumberParticles(L, s);

    // Vector used (and labelled w.r.t particles location)
    vector<double> pos_array, type_arr;
    vector<vector<int>> cell_matrix(Nx * Ny * Nz);

    // Variables defined to used "export.cpp"
    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;

    /*---------------------------- SPH ALGORITHM  ----------------------------*/

    // Initialise variables of the problem
    initializeRho(params, pos, rho);
    initializeMass(params, rho, mass);
    initializeVelocity(params, u);
    initializeViscosity(params, artificial_visc_matrix);
    setPressure(params, p, rho); 
    setSpeedOfSound(params, c, rho);

    // Check if the chose, timeStep is small enough
    checkTimeStep(params, pos, c, neighbours_matrix, artificial_visc_matrix);
    

    for (int t = 0; t < params.nstepT; t++){

        // Apply the linked-list algorithm
        sorted_list(params, cell_matrix, neighbours_matrix, gradW_matrix, pos); 

        // Compute ∇_a(W_ab) for all particles
        gradW(gradW_matrix, neighbours_matrix, 
              pos_array,
              nb_moving_part, h, Nx, Ny, Nz,
              PRINT); 

        // Compute pressure for all particles
        setPressure(p_array, rho_array, 
                    nb_moving_part, rho_0, c_0, R, T, M, gamma, 
                    state_equation_chosen, PRINT); 

        // After updates, need to check if timeStep is still small enough
        checkTimeStep(params, pos, c, neighbours_matrix, artificial_visc_matrix);


        // Clear matrices and reset arrays to 0
        clearAllVectors(params, artificial_visc_matrix, neighbours_matrix, cell_matrix, gradW_matrix,
                        drhodt, dudt);

        // Compute artificial viscosity Π_ab for all particles
        setArtificialViscosity(artificial_visc_matrix, neighbours_matrix,
                               c_array, pos_array, rho_array, u_array, 
                               nb_moving_part, t, alpha, beta, rho_0, c_0, gamma, R, T, M, h,
                               state_equation_chosen, PRINT); 

        // Compute D(rho)/Dt for all particles
        continuityEquation(neighbours_matrix, gradW_matrix, 
                           pos_array, u_array, drhodt_array, rho_array, mass_array,
                           nb_moving_part, h,
                           PRINT); 

        // Compute D(u)/Dt for all particles
        momentumEquation(neighbours_matrix, gradW_matrix, artificial_visc_matrix,
                         mass_array, dudt_array, rho_array, p_array,
                         nb_moving_part,
                         rho_0, c_0, gamma, R, T, M, g, 
                         state_equation_chosen, PRINT); 

        // Update density, velocity and position for each particle (Euler explicit scheme)
        for (size_t pos = 0; pos < nb_tot_part; pos++){

            rho_array[pos] += dt * drhodt_array[pos];

            for (size_t cord = 0; cord < 3; cord++)
            {

                pos_array[3 * pos + cord] += dt * u_array[3 * pos + cord];
                u_array[3 * pos + cord] += dt * dudt_array[3 * pos + cord];
            }
        }
        for(size_t i = 0 ; i < gradW_matrix.size(); i++){
            for(size_t j = 0 ; j < gradW_matrix[i].size(); j ++){
                grad_sum[i] += abs(gradW_matrix[i][j]);
            }
        }
        
        clearAllVectors(artificial_visc_matrix, neighbours_matrix,
                        cell_matrix, gradW_matrix,
                        PRINT);

        if(t % nsave == 0){
            export_particles("../../output/sph", t, pos_array, scalars, vectors);
        }
        for(size_t i = 0 ; i<drhodt_array.size(); i ++ )
        {
            drhodt_array[i] = 0.0;
        }
        for(size_t i = 0 ; i<dudt_array.size(); i ++ )
        {
            dudt_array[i] = 0.0;
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    std::cout << "duration: " << delta_t << "s.\n";

    std::cout << "\nSimulation done." << std::endl;
    return 0;
}
