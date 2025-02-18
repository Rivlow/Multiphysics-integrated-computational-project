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
#include "NavierStokes.h"
#include "initialize.h"
#include "export.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"
#include "time_integration.h"
#include "data_store.h"
#include "surface_tension.h"


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

    //omp_set_num_threads(12);
   
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
    string state_equation;
    string state_initial_condition;
    string scheme_integration;
    string kernel;
    string name_file = data["name_file"];

    GeomData geomParams;
    SimulationData simParams;
    ThermoData thermoParams;

    if (data["omp"]["chose_nb_of_threads"])
        omp_set_num_threads(data["omp"]["nb_of_threads"]);

    // Select desired simulation inputs
    getKey(data, state_equation, state_initial_condition, 
           scheme_integration, kernel);

    // Structure to store parameters
    initializeStruct(data, state_equation, state_initial_condition, 
                      scheme_integration, kernel, geomParams, simParams, thermoParams);

    // Make sure output folders are present and empty for current simulation
    createOutputFolder();
    clearOutputFiles();

    /*------------------- INITIALIZATION OF VARIABLES USED -------------------*/    

    // Initialize particles
    int MP_count = 0, FP_count = 0, GP_count = 0;
    vector<double> pos;
    vector<double> type;
    
    meshCube(geomParams, simParams, pos, type, MP_count, FP_count); 
    meshPostProcess(geomParams, simParams, pos, type, GP_count);
    if (!checkParticleGeneration(pos, simParams))
        exit(1);
    
    int nb_tot_part = pos.size()/3;
    vector<double> mass(simParams.nb_tot_part + GP_count, 0),
                   u(3 * (simParams.nb_tot_part + GP_count), 0),
                   drhodt(simParams.nb_tot_part + GP_count, 0), 
                   rho(simParams.nb_tot_part + GP_count, 0),
                   dudt(3 * (simParams.nb_tot_part + GP_count), 0), 
                   p(simParams.nb_tot_part + GP_count, 0),
                   c(simParams.nb_tot_part + GP_count, 0), 
                   grad_sum(simParams.nb_tot_part + GP_count, 0),
                   track_particle(simParams.nb_tot_part + GP_count, 0),
                   nb_neighbours(simParams.nb_tot_part + GP_count, 0),
                   N(simParams.nb_tot_part + GP_count, 0),
                   colour(simParams.nb_tot_part + GP_count, 0),
                   normal(3*(simParams.nb_tot_part + GP_count), 0),
                   acc_vol(3*(simParams.nb_tot_part + GP_count), 0),
                   R(simParams.nb_tot_part + GP_count, 0),
                   L(simParams.nb_tot_part + GP_count, 0),
                   Kappa(simParams.nb_tot_part + GP_count, 0),
                   dot_product(simParams.nb_tot_part + GP_count, 0);

    vector<vector<int>> cell_matrix(geomParams.Nx * geomParams.Ny * geomParams.Nz);
    vector<vector<double>> gradW(simParams.nb_tot_part + GP_count),
                           W(simParams.nb_tot_part + GP_count),
                           gradW_st(simParams.nb_tot_part + GP_count),
                           W_st(simParams.nb_tot_part + GP_count),
                           viscosity(simParams.nb_tot_part + GP_count);

    vector<int> neighbours(100*(simParams.nb_tot_part + GP_count));

    // Variables defined to used "export.cpp"
    map<string, vector<double> *> scalars;
    map<string, vector<double> *> vectors;
    scalars["type"] = &type;
    scalars["track_particle"] = &track_particle;
    scalars["nb_neighbours"] = &nb_neighbours;
    scalars["mass"] = &mass;
    scalars["rho"] = &rho;
    scalars["p"] = &p;
    scalars["drhodt"] = &drhodt;
    vectors["position"] = &pos;
    vectors["u"] = &u;
    vectors["dudt"] = &dudt;
    scalars["colour"] = &colour;
    scalars["R"] = &R;
    scalars["L"] = &L;
    scalars["N"] = &N;
    scalars["track_particle"] = &track_particle;
    scalars["Kappa"] = &Kappa;
    vectors["normal"]= &normal;
    vectors["acc_vol"]= &acc_vol;

    printParams(data, geomParams, thermoParams, simParams,
                state_equation, state_initial_condition, scheme_integration,
                MP_count, FP_count, GP_count, nb_tot_part);

    /*---------------------------- SPH ALGORITHM  ----------------------------*/

    // Initialize variables of the problem
    initRho(thermoParams, simParams, pos, rho);
    initMass(geomParams, simParams, rho, mass);
    initVelocity(thermoParams, simParams, u);
    setPressure(geomParams, thermoParams, simParams, p, rho); 
    setSpeedOfSound(geomParams, thermoParams, simParams, c, rho);
    initKernelCoef(geomParams, simParams);


    // If needed, compare linked-list and naive searching neighbour algorithm
    if (simParams.comparison_algorithm){
        compareAlgo(geomParams, simParams, cell_matrix, neighbours, gradW,
                    W, viscosity, nb_neighbours, type, pos);
        return 0;
    }

    auto t_mid = chrono::high_resolution_clock::now();
    double sim_time = 0.0;
    int iteration = 0;
    vector<double> vec_time(simParams.nstepT/simParams.nsave, 0.0);

    for (int t = 0; t < simParams.nstepT; t++){

        simParams.t = t;

        // Apply the linked-list algorithm
        sortedList(geomParams, simParams, cell_matrix, neighbours,
                   gradW, W, viscosity, nb_neighbours, type, pos);
    
        // Compute ∇_a(W_ab) for all particles (moving and fixed)
        computeGradW(geomParams, simParams, gradW, W, neighbours, nb_neighbours, pos);

        // Update density, velocity and position (Euler explicit or RK22 scheme)
        updateVariables(geomParams, thermoParams, simParams, pos, u, rho, drhodt, c, p, dudt, mass, 
                        viscosity, gradW, W, neighbours, nb_neighbours, type, colour, R, L, N, normal, acc_vol, track_particle, Kappa, dot_product);

        sim_time += simParams.dt;

        // Save data each "nsave" iterations
        if(t % simParams.nsave == 0){
    
            if (geomParams.post_process_do)
                extractData(geomParams, simParams, thermoParams, scheme_integration, name_file, pos, p, mass, u, neighbours, nb_neighbours, rho);
            if (geomParams.following_part_bool || geomParams.following_part_max || geomParams.following_part_min)
                follow_part_data(geomParams,simParams, name_file, p, rho, pos, u);
               
            
            export_particles("../../output/sph", t, pos, scalars, vectors, false);
            auto t_act = chrono::high_resolution_clock::now();
            double elapsed_time = double(chrono::duration_cast<chrono::duration<double>>(t_act - t_mid).count());
            progressBar(double(t)/double(simParams.nstepT), elapsed_time); 
            vec_time[iteration++] = sim_time; 
        }

        // Clear matrices and reset arrays to 0
        clearAllVectors(simParams, viscosity, neighbours,
                        cell_matrix, gradW, W, drhodt, dudt, colour, R, N, normal, acc_vol, track_particle, Kappa, dot_product);

    }

    writing_time(vec_time, name_file,simParams);
    auto t1 = chrono::high_resolution_clock::now();
    auto delta_t = chrono::duration_cast<chrono::duration<double>>(t1 - t0).count();
    cout << "\nDuration: " << delta_t << "s (" << int((MP_count+FP_count)/delta_t)*simParams.nstepT<<" [Particles*Updates/s]).\n";
    cout << "Simulation done." << endl;

    return 0;
}
