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

        if (argc != 2){

            cout << "\nusage: " << argv[0] << " <param.json>\n\n";
            return EXIT_FAILURE;
        }


    /*---------------------- INPUT PARAMETERS FROM JSON FILES --------------------------*/

    ifstream inputf(argv[1]);
    json data = json::parse(inputf);

    //cout << argv[1] << ":\n"
    //          << data.dump(4) << endl; // print input data to screen

    createOutputFolder();
    clearOutputFiles();

    string state_equation;
    string state_initial_condition;
    string schemeIntegration;

    // Structure to store parameters
    getKey(data, state_equation, state_initial_condition, 
           schemeIntegration);

    cout << "kappa" << data["kappa"] << endl;
    cout << "s" << data["s"] << endl;
    cout << "o" << data["domain"]["o"] << endl;
    cout << "L" << data["domain"]["L"] << endl;
    cout << "o_d" << data["domain"]["o_d"] << endl;
    cout << "L_d" << data["domain"]["L_d"]<< endl;
    cout << "do" << data["post_process"]["do"] << endl;
    cout << "xyz_init" << data["post_process"]["xyz_init"] << endl;
    cout << "xyz_end" << data["post_process"]["xyz_end"] << endl;
    cout << "matrix_long" << data["domain"]["matrix_long"] << endl;
    cout << "matrix_orig" << data["domain"]["matrix_orig"] << endl;
    cout << "vector_type" << data["domain"]["vector_type"] << endl;


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
    cout << "GeomData initialized" << endl;

    ThermoData thermoParams = {
        data["thermo"]["c_0"],
        data["thermo"]["rho_moving"],
        data["thermo"]["rho_fixed"],
        data["thermo"]["rho_0"],
        data["thermo"]["M"],
        data["thermo"]["T"],
        data["thermo"]["gamma"],
        data["thermo"]["R"],
    };
    cout << "ThermoData initialized" << endl;


    SimulationData simParams = {

        data["simulation"]["nstepT"],
        data["simulation"]["nsave"],
        data["simulation"]["dt"],
        data["simulation"]["theta"],
        data["simulation"]["alpha"],
        data["simulation"]["beta"],
        data["simulation"]["alpha_st"],
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
    vector<vector<int>> cell_matrix(geomParams.Nx * geomParams.Ny * geomParams.Nz);

    // Initialize particles
    int MP_count = 0, FP_count = 0, GP_count = 0;
    meshcube(geomParams, simParams, pos, type, MP_count, FP_count); 
    meshPostProcess(geomParams, simParams, pos, type, GP_count);
    
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

    cout << "#===============================#" << endl;
    cout << "# General simulation parameters #" << endl;
    cout << "#===============================#" << "\n" << endl;
    cout << "s =" << geomParams.s << endl;
    cout << "kappa = " << geomParams.kappa << endl;
    cout << "h = " << geomParams.h << endl;
    cout << "nstepT = " << simParams.nstepT << endl;
    cout << "nsave = " << simParams.nsave << endl;
    cout << "dt = " << simParams.dt << endl;
    cout << "theta = " << simParams.theta <<endl;
    cout << "alpha (artificial viscosity): " << simParams.alpha << endl;
    cout << "alpha (surface tension): " << simParams.alpha_st << endl;
    cout << "beta (artificial visocity): " << simParams.beta  << endl;
    cout << "state equation : " << state_equation << endl;
    cout << "state initial condition : " << state_initial_condition << "\n" << endl;

    cout << "#==================#" << endl;
    cout << "# Domain variables #" << endl;
    cout << "#==================#" << "\n" << endl;
    cout << "Radius of neighbourhood (kappa * h) = " << geomParams.kappa * geomParams.h << endl;
    cout << "Number of cells in each direction (Nx, Ny, Nz) = (" << geomParams.Nx << ", " << geomParams.Ny << ", " << geomParams.Nz << ")" << endl;
    cout << "Number of fluid particles = " << MP_count << endl;
    cout << "Number of fixed particles = " << FP_count << endl;
    cout << "Number of post process particles = " << GP_count << endl;
    cout << "Total number of particles = " << nb_tot_part << "\n" << endl;

    cout << "#==================#" << endl;
    cout << "# Thermo variables #" << endl;
    cout << "#==================#" << "\n" << endl; 
    cout << "Initial speed of sound (c_0) = " << thermoParams.c_0 << endl;
    cout << "Moving particle density (rho_moving) = " << thermoParams.rho_moving << endl;
    cout << "Fixed particle density (rho_fixed) = " << thermoParams.rho_fixed << endl;
    cout << "Initial density (rho_0) = " << thermoParams.rho_0 << endl;
    cout << "Molar mass (M) = " << thermoParams.M << endl;
    cout << "Heat capacity ration (gamma) = " << thermoParams.gamma << endl;
    cout << "Ideal gaz constant (R) = " << thermoParams.R << "\n" << endl;

    
    

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

        // Apply the linked-list algorithm
        sortedList(geomParams, simParams, cell_matrix, neighbours, 
                   gradW_matrix, pi_matrix, nb_neighbours, type, pos); 

        // Compute ∇_a(W_ab) for all particles
        gradW(geomParams, simParams, gradW_matrix, neighbours, nb_neighbours, pos); 

        // Update density, velocity and position (Euler explicit or RK22 scheme)
        updateVariables(geomParams, thermoParams, simParams, pos, u, rho, drhodt, c, p, dudt, mass, 
                        pi_matrix, gradW_matrix, neighbours, nb_neighbours);

        // Check if timeStep is small enough
        checkTimeStep(geomParams, thermoParams, simParams, pos, u, c,
                      neighbours, nb_neighbours, pi_matrix);

        // Save data each "nsave" iterations
        if(t % simParams.nsave == 0){
                if (geomParams.post_process_do)
                    extractData(geomParams, simParams, thermoParams, pos, p, mass, neighbours, nb_neighbours);
                
            export_particles("../../output/sph", t, pos, scalars, vectors, false);

            auto t_act = chrono::high_resolution_clock::now();
            double elapsed_time = double(chrono::duration_cast<chrono::duration<double>>(t_act - t0).count());
            progressBar(double(t)/double(simParams.nstepT), elapsed_time);

                        
            // Calculer le temps écoulé depuis le début de la simulation
            
            
        }


        // Clear matrices and reset arrays to 0
        clearAllVectors(simParams, pi_matrix, neighbours,
                        cell_matrix, gradW_matrix, drhodt, dudt);

    }

    auto t1 = chrono::high_resolution_clock::now();
    auto delta_t = chrono::duration_cast<chrono::duration<double>>(t1 - t0).count();
    cout << "\nduration: " << delta_t << "s (" << int((MP_count+FP_count)/delta_t)<<" particles/s).\n";
    cout << "Simulation done." << endl;

    return 0;
}
