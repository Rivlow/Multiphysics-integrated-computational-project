#include <stdio.h>
#include <string.h>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>

#include "generate_particle.h"
#include "find_neighbours.h"
#include "gradient.h"
#include "initialize.h"
#include "export.h"
#include "Kernel_functions.h"
#include "tools.h"
#include "structure.h"


#ifdef _OPENMP
#include <omp.h>
#endif
#include <chrono>

#include "nlohmann/json.hpp"
using json = nlohmann::json;
using namespace std;


int main(int argc, char *argv[])
{

    /*---------------- SETTING COMPILATION PARAMETERS/FOLDER NEEDED --------------------*/

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


    /*---------------------- INPUT PARAMETERS FROM JSON FILES --------------------------*/

    std::ifstream inputf(argv[1]);
    json data = json::parse(inputf);

    std::cout << argv[1] << ":\n"
              << data.dump(4) << std::endl; // print input data to screen


    std::string state_equation;
    std::string state_initial_condition;

    for (auto &it : data["stateEquation"].items())
    {
        if (it.value() == true)
        {
            state_equation = it.key();
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


    // Structure to store all parameters used later 
    const SimulationData params = {

        data["kappa"],
        data["nstepT"],
        data["nsave"],
        data["dt"],
        data["s"],
        1.2 * params.s,
        data["o"],
        data["L"],
        data["o_d"],
        data["L_d"],
        data["u_init"],
        int(params.L_d[0] / (params.kappa * params.h)),
        int(params.L_d[1] / (params.kappa * params.h)),
        int(params.L_d[2] / (params.kappa * params.h)),

        data["alpha"],
        data["beta"],

        data["c_0"],
        data["rho_moving"],
        data["rho_fixed"],
        data["rho_0"],
        data["M"],
        data["T"],
        data["gamma"],
        8.314, // [J/(K.mol)]
        -9.81, // [m/s²]
        state_equation,
        state_initial_condition,
        data["print_debug"],
        evaluateNumberParticles(params),
    };


    /*------------------- INITIALIZATION OF VARIABLES USED -------------------*/

    // Vector used (and labelled w.r.t particles location)
    vector<double> pos_array, type_arr;
    vector<vector<int>> cell_matrix(params.Nx * params.Ny * params.Nz);

    // Initialization of the particles (moving and fixed)
    meshcube(params, pos_array, type_arr);
    meshBoundary(params, pos_array, type_arr);
    int nb_tot_part = pos_array.size()/3;

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
    std::vector<double> nvoisins(nb_tot_part, 0.0); 

    // Variables defined to used "export.cpp"
    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;


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



    /*---------------------------- SPH ALGORITHM  ----------------------------*/

    cout << "state equation chosen : " << state_equation << " \n" << endl;
    cout << "kappa * h =" << params.kappa * params.h << endl;
    cout << "(Nx, Ny, Nz) = (" << params.Nx << ", " << params.Ny << ", " << params.Nz << ")" << std::endl;
    cout << "b_moving_part = " << params.nb_moving_part << std::endl;
    cout << "nb_tot_part = " << nb_tot_part << std::endl;
    cout << "s=" << params.s << std::endl;
    cout << "kappa=" << params.kappa << std::endl;
    cout << "h=" << params.h << std::endl;

    

    initializeRho(params, pos_array, rho_array);
    initializeMass(params, rho_array, mass_array);
    initializeVelocity(params, u_array);
    initializeViscosity(params, artificial_visc_matrix);


    for (int t = 0; t < params.nstepT; t++)
    {

        // compute dt_max = min_a h_a/F
        // double dt_max = sqrt(kappa * h / abs(g)); // CFL condition
        // std::cout << "dt_max = " << dt_max << "    dt = " << dt << std::endl;

        // Apply the linked-list algorithm
        sorted_list(params, cell_matrix, neighbours_matrix, pos_array); 


        // Compute ∇_a(W_ab) for all particles
        gradW(params, gradW_matrix, neighbours_matrix, 
              pos_array); 


        // Compute D(rho)/Dt for all particles
        continuityEquation(params, neighbours_matrix, gradW_matrix, 
                           pos_array, u_array, drhodt_array, rho_array, mass_array); 


        // Compute D(u)/Dt for all particles
        momentumEquation(params, t, neighbours_matrix, gradW_matrix, artificial_visc_matrix,
                         mass_array, dudt_array, rho_array, p_array, c_array, pos_array, u_array); 
        

        // Update density, velocity and position for each particle (Euler explicit scheme)
        update(params, pos_array, u_array, rho_array, drhodt_array, dudt_array);
        
        // Clear matrices and reset arrays to 0
        clearAllVectors(params, artificial_visc_matrix, neighbours_matrix, cell_matrix, gradW_matrix,
                        drhodt_array, dudt_array);


        if(t % params.nsave == 0){
            export_particles("../../output/sph", t, pos_array, scalars, vectors);
        }

    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    std::cout << "duration: " << delta_t << "s.\n";

    std::cout << "\n Simulation done." << std::endl;
    return 0;
}
