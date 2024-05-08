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

    auto t0 = chrono::high_resolution_clock::now();

    #ifdef _OPENMP
        cout << "OpenMP available: OMP_NUM_THREADS=" << omp_get_max_threads() << "\n";
    #else
        cout << "OpenMP not available.\n";
    #endif

    #ifdef NDEBUG
        // code has been configured with "cmake -DCMAKE_BUILD_TYPE=Release .."
        cout << "code built in RELEASE mode.\n";
    #else
        // code has been configured with "cmake .."
        cout << "code built in DEBUG mode.\n";
    #endif

        if (argc != 2)
        {
            cout << "\nusage: " << argv[0] << " <param.json>\n\n";
            return EXIT_FAILURE;
        }


    /*---------------------- INPUT PARAMETERS FROM JSON FILES --------------------------*/

    ifstream inputf(argv[1]);
    json data = json::parse(inputf);

    cout << argv[1] << ":\n"
              << data.dump(4) << endl; // print input data to screen

    createOutputFolder();
    clearOutputFiles();

    string state_equation;
    string state_initial_condition;
    string schemeIntegration;

    // Structure to store parameters
    getKey(data, state_equation, state_initial_condition, 
           schemeIntegration);

    GeomData geomParams = {

        data["simulation"]["kappa"],
        data["simulation"]["s"],
        1.2 * geomParams.s,
        data["domain"]["o_d"],
        data["domain"]["L_d"],
        data["domain"]["matrix_long"],
        data["domain"]["matrix_orig"],
        data["domain"]["vector_type"],
        data["post_process"]["xyz_init"],
        data["post_process"]["xyz_end"],
        data["post_process"]["do"],
        int(geomParams.L_d[0] / (geomParams.kappa * geomParams.h)),
        int(geomParams.L_d[1] / (geomParams.kappa * geomParams.h)),
        int(geomParams.L_d[2] / (geomParams.kappa * geomParams.h)),

    };

    ThermoData thermoParams = {

        data["simulation"]["alpha"],
        data["simulation"]["beta"],
        data["thermo"]["c_0"],
        data["thermo"]["rho_moving"],
        data["thermo"]["rho_fixed"],
        data["thermo"]["rho_0"],
        data["thermo"]["M"],
        data["thermo"]["T"],
        data["thermo"]["gamma"],
        data["thermo"]["R"],
    };

    SimulationData simParams = {

        data["simulation"]["nstepT"],
        data["simulation"]["nsave"],
        data["simulation"]["dt"],
        data["simulation"]["theta"],
        schemeIntegration,
        data["thermo"]["u_init"],
        state_equation,
        state_initial_condition,
        false,
        true,
        data["condition"]["print_debug"],
        evaluateNumberParticles(geomParams),
        0,
        0,
    };




    /*------------------- INITIALIZATION OF VARIABLES USED -------------------*/

    // Vector used (and labelled w.r.t particles location)
    vector<double> pos, type;
    int Nx = geomParams.Nx, Ny = geomParams.Ny, Nz = geomParams.Nz; 
    vector<vector<int>> cell_matrix(Nx * Ny * Nz);

    // Initialize particles
    meshcube(geomParams, simParams, pos, type); 
    if (geomParams.post_process_do) 
        meshPostProcess(geomParams, simParams, pos, type);
    int nb_tot_part = pos.size()/3;

    vector<double> mass(nb_tot_part, 0), u(3 * nb_tot_part, 0),
                   drhodt(nb_tot_part, 0), rho(nb_tot_part, 0),
                   dudt(3 * nb_tot_part, 0), p(nb_tot_part, 0),
                   c(nb_tot_part, 0), grad_sum(nb_tot_part, 0);

    vector<vector<double>> pi_matrix(nb_tot_part), gradW_matrix(nb_tot_part);
    vector<int> neighbours(100*nb_tot_part);
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
    cout << "kappa * h = " << geomParams.kappa * geomParams.h << endl;
    cout << "(Nx, Ny, Nz) = (" << geomParams.Nx << ", " << geomParams.Ny << ", " << geomParams.Nz << ")" << endl;
    cout << "nb_moving_part = " << simParams.nb_moving_part << endl;
    cout << "nb_tot_part = " << nb_tot_part << endl;
    cout << "s =" << geomParams.s << endl;
    cout << "kappa = " << geomParams.kappa << endl;
    cout << "h = " << geomParams.h << endl;

    /*---------------------------- SPH ALGORITHM  ----------------------------*/

    // Initialize variables of the problem
    initRho(thermoParams, simParams, pos, rho);
    initMass(geomParams, simParams, rho, mass);
    initVelocity(thermoParams, simParams, u);
    initViscosity(simParams, pi_matrix);
    setPressure(geomParams, thermoParams, simParams, p, rho); 
    setSpeedOfSound(geomParams, thermoParams, simParams, c, rho);
    

    for (int t = 0; t < simParams.nstepT; t++){

        simParams.t = t;

        // Check if timeStep is small enough
        checkTimeStep(geomParams, thermoParams, simParams, pos, u, c,
                      neighbours, nb_neighbours, pi_matrix);

        // Apply the linked-list algorithm
        sortedList(geomParams, simParams, cell_matrix, neighbours, 
                   gradW_matrix, pi_matrix, nb_neighbours, type, pos); 

        // Compute ∇_a(W_ab) for all particles
        gradW(geomParams, simParams, gradW_matrix, neighbours, nb_neighbours, pos); 

        // Update density, velocity and position (Euler explicit or RK22 scheme)
        updateVariables(geomParams, thermoParams, simParams, pos, u, rho, drhodt, c, p, dudt, mass, 
                        pi_matrix, gradW_matrix, neighbours, nb_neighbours);

        // Save data each "nsave" iterations
        if(t % simParams.nsave == 0){
                if (geomParams.post_process_do)
                    extractData(geomParams, simParams, thermoParams, pos, p, mass, neighbours, nb_neighbours);
                
            export_particles("../../output/sph", t, pos, scalars, vectors, false);
            progresssBar(double(t)/double(simParams.nstepT) * 100, t, simParams.nstepT);
        }


        // Clear matrices and reset arrays to 0
        clearAllVectors(simParams, pi_matrix, neighbours,
                        cell_matrix, gradW_matrix, drhodt, dudt);

    }

    auto t1 = chrono::high_resolution_clock::now();
    auto delta_t = chrono::duration_cast<chrono::duration<double>>(t1 - t0).count();
    cout << " \n duration: " << delta_t << "s.\n";
    cout << "Simulation done." << endl;

    return 0;
}
