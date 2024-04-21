#include <stdio.h>
#include <string.h>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <chrono>

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

    ifstream inputf(argv[1]);
    json data = json::parse(inputf);

    std::cout << argv[1] << ":\n"
              << data.dump(4) << std::endl; // print input data to screen

    createOutputFolder();
    clearOutputFiles();

    string state_equation;
    string state_initial_condition;
    string schemeIntegration;
    vector<string> walls_chose;

    // Structure to store parameters
    getKey(data, state_equation, state_initial_condition, 
           schemeIntegration, walls_chose);


    GeomData geomParams = {

        data["kappa"],
        data["s"],
        1.2 * geomParams.s,
        data["o"],
        data["L"],
        data["o_d"],
        data["L_d"],
        int(geomParams.L_d[0] / (geomParams.kappa * geomParams.h)),
        int(geomParams.L_d[1] / (geomParams.kappa * geomParams.h)),
        int(geomParams.L_d[2] / (geomParams.kappa * geomParams.h)),

        data["matrixLong"],
        data["matrixOrig"],
        data["vectorType"],
    };

    ThermoData thermoParams = {

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
    };

    SimulationData simParams = {

        data["nstepT"],
        data["nsave"],
        data["dt"],
        data["theta"],
        schemeIntegration,
        data["data_store"]["name"],
        data["data_store"]["init"],
        data["data_store"]["end"],
        data["data_store"]["do"],
        data["u_init"],
        state_equation,
        state_initial_condition,
        data["print_debug"],
        evaluateNumberParticles(geomParams),
    };


    /*------------------- INITIALIZATION OF VARIABLES USED -------------------*/

    // Vector used (and labelled w.r.t particles location)
    vector<double> pos, type;

    int Nx = geomParams.Nx, Ny = geomParams.Ny, Nz = geomParams.Nz; 

    vector<vector<int>> cell_matrix(Nx * Ny * Nz);

    // Initialization of the particles
    meshcube(geomParams, pos, type); // moving
    //meshBoundary(geomParams, pos, type); // fixed
    int nb_tot_part = pos.size()/3;

    vector<double> mass(nb_tot_part), u(3 * nb_tot_part),
                   drhodt(nb_tot_part), rho(nb_tot_part),
                   dudt(3 * nb_tot_part, 0.0), p(nb_tot_part),
                   c(nb_tot_part), grad_sum(nb_tot_part);

    vector<vector<double>> pi_matrix(nb_tot_part), gradW_matrix(nb_tot_part);
    vector<vector<int>> neighbours_matrix(nb_tot_part, vector<int>(100));
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
    cout << "kappa * h =" << geomParams.kappa * geomParams.h << endl;
    cout << "(Nx, Ny, Nz) = (" << geomParams.Nx << ", " << geomParams.Ny << ", " << geomParams.Nz << ")" << std::endl;
    cout << "b_moving_part = " << simParams.nb_moving_part << std::endl;
    cout << "nb_tot_part = " << nb_tot_part << std::endl;
    cout << "s=" << geomParams.s << std::endl;
    cout << "kappa=" << geomParams.kappa << std::endl;
    cout << "h=" << geomParams.h << std::endl;


    //printArray(type, type.size(),"type");

    /*---------------------------- SPH ALGORITHM  ----------------------------*/

    // Initialise variables of the problem
    initializeRho(thermoParams, simParams, pos, rho);
    initializeMass(geomParams, simParams, rho, mass);
    initializeVelocity(thermoParams, simParams, u);
    initializeViscosity(simParams, pi_matrix);
    setPressure(geomParams, thermoParams, simParams, p, rho); 
    setSpeedOfSound(geomParams, thermoParams, simParams, c, rho);

    for (int t = 0; t < simParams.nstepT; t++){
        
        // Check if timeStep is small enough
        checkTimeStep(geomParams, thermoParams, simParams, t, pos, c, neighbours_matrix, nb_neighbours, pi_matrix);

        // Apply the linked-list algorithm
        sortedList(geomParams, simParams, cell_matrix, neighbours_matrix, gradW_matrix, 
                    pi_matrix, nb_neighbours, pos); 
        //printMatrix(neighbours_matrix, neighbours_matrix.size(), "neig");
        // Compute ∇_a(W_ab) for all particles
        gradW(geomParams, simParams, gradW_matrix, neighbours_matrix, nb_neighbours, pos); 

        // Update density, velocity and position (Euler explicit or RK22 scheme)
        updateVariables(geomParams, thermoParams, simParams, t, pos, u, rho, drhodt, c, p, dudt, mass, 
                        pi_matrix, gradW_matrix, neighbours_matrix, nb_neighbours);
        //printArray(dudt,3*simParams.nb_moving_part,"dudt");
        // Save data each "nsave" iterations
        if(t % simParams.nsave == 0){

            if (simParams.data_do){extractData(simParams, pos, u, dudt, rho, drhodt, c, p, mass);}
            export_particles("../../output/sph", t, pos, scalars, vectors);

        }

        // Clear matrices and reset arrays to 0
        clearAllVectors(simParams, pi_matrix, neighbours_matrix,
                        cell_matrix, gradW_matrix, drhodt, dudt);

         //printArray(dudt,3*simParams.nb_moving_part,"dudt_clear");
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto delta_t = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    std::cout << "duration: " << delta_t << "s.\n";

    std::cout << "\n Simulation done." << std::endl;
    return 0;
}
