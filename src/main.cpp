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
#include "Kernel.h"
#include "tools.h"
#include "structure.h"
#include "integration.h"
#include "data_store.h"



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

    DomainParams domainParams = {
        data["domain"]["shape"],
        {},
        data["domain"]["particle_layers"],

    };

    for (auto& wall : data["domain"]["walls_used"].items()) {
        if (wall.value().get<bool>()) {
            domainParams.walls_used.push_back(wall.key());
        }
    }

     SimulationData params = {

        data["kappa"],
        data["nstepT"],
        data["nsave"],
        data["dt"],
        data["theta"],
        data["schemeIntegration"],
        data["data_store"]["name"],
        data["data_store"]["init"],
        data["data_store"]["end"],
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
        domainParams,

    };


    /*------------------- INITIALIZATION OF VARIABLES USED -------------------*/

    // Vector used (and labelled w.r.t particles location)
    vector<double> pos, type;
    vector<vector<int>> cell_matrix(params.Nx * params.Ny * params.Nz);

    // Initialization of the particles
    meshcube(params, pos, type); // moving
    meshBoundary(params, pos, type); // fixed
    int nb_tot_part = pos.size()/3;

    vector<double> mass(nb_tot_part), u(3 * nb_tot_part),
                   drhodt(nb_tot_part), rho(nb_tot_part),
                   dudt(3 * nb_tot_part, 0.0), p(nb_tot_part),
                   c(nb_tot_part), grad_sum(nb_tot_part);

    vector<vector<double>> pi_matrix(nb_tot_part), gradW_matrix(nb_tot_part);
    vector<vector<int>> neighbours_matrix(nb_tot_part);
    vector<double> nb_neighbours(nb_tot_part, 0.0); 

    // Variables defined to used "export.cpp"
    map<string, vector<double> *> scalars;
    map<string, vector<double> *> vectors;


    scalars["type"] = &type;
    scalars["mass"] = &mass;
    scalars["rho"] = &rho;
    scalars["p"] = &p;
    scalars["drhodt"] = &drhodt;
    scalars["nb_neighbours"] = &nb_neighbours;
    scalars["grad_sum"] = &grad_sum;
    vectors["position"] = &pos;
    vectors["u"] = &u;
    vectors["dudt"] = &dudt;

    cout << "state equation chosen : " << state_equation << " \n" << endl;
    cout << "kappa * h =" << params.kappa * params.h << endl;
    cout << "(Nx, Ny, Nz) = (" << params.Nx << ", " << params.Ny << ", " << params.Nz << ")" << std::endl;
    cout << "b_moving_part = " << params.nb_moving_part << std::endl;
    cout << "nb_tot_part = " << nb_tot_part << std::endl;
    cout << "s=" << params.s << std::endl;
    cout << "kappa=" << params.kappa << std::endl;
    cout << "h=" << params.h << std::endl;




    /*---------------------------- SPH ALGORITHM  ----------------------------*/

    // Initialise variables of the problem
    initializeRho(params, pos, rho);
    initializeMass(params, rho, mass);
    initializeVelocity(params, u);
    initializeViscosity(params, pi_matrix);
    setPressure(params, p, rho); 
    setSpeedOfSound(params, c, rho);

    for (int t = 0; t < params.nstepT; t++){

        // Check if timeStep is small enough
        checkTimeStep(params, t, pos, c, neighbours_matrix, pi_matrix);

        // Apply the linked-list algorithm
        sorted_list(params, cell_matrix, neighbours_matrix, gradW_matrix, 
                    pi_matrix, nb_neighbours, pos); 

        // Compute ∇_a(W_ab) for all particles
        gradW(params, gradW_matrix, neighbours_matrix, pos); 

        // Update density, velocity and position (Euler explicit or RK22 scheme)
        updateVariables(params, t, pos, u, rho, drhodt, c, p, dudt, mass, 
                        pi_matrix, gradW_matrix, neighbours_matrix);

        if(t % params.nsave == 0){
            export_particles("../../output/sph", t, pos, scalars, vectors);
            extractData(params, pos, u, dudt, rho, drhodt, c, p, mass);

        }

        // Clear matrices and reset arrays to 0
        clearAllVectors(params, pi_matrix, neighbours_matrix,
                        cell_matrix, gradW_matrix, drhodt, dudt);


    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    std::cout << "duration: " << delta_t << "s.\n";

    std::cout << "\n Simulation done." << std::endl;
    return 0;
}
