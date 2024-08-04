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

    cout << "kappa = " << data["simulation"]["kappa"] << endl;
    cout << "s = " << data["simulation"]["s"] << endl;
    cout << "o_d = " << data["domain"]["o_d"] << endl;
    cout << "L_d = " << data["domain"]["L_d"]<< endl;
    cout << "do = " << data["post_process"]["do"] << endl;
    cout << "xyz_init = " << data["post_process"]["xyz_init"] << endl;
    cout << "xyz_end = " << data["post_process"]["xyz_end"] << endl;
    cout << "matrix_long = " << data["domain"]["matrix_long"] << endl;
    cout << "matrix_orig = " << data["domain"]["matrix_orig"] << endl;
    cout << "vector_type = " << data["domain"]["vector_type"] << endl;


    GeomData geomParams = {

        data["simulation"]["kappa"],
        data["simulation"]["s"],
        1.2 * geomParams.s,
        data["domain"]["o_d"],
        data["domain"]["L_d"],
        data["domain"]["matrix_long"],
        data["domain"]["matrix_orig"],
        data["domain"]["vector_type"],
        data["domain"]["sphere"]["do"],
        data["domain"]["sphere"]["radius"],
        data["post_process"]["xyz_init"],
        data["post_process"]["xyz_end"],
        data["post_process"]["do"],
        int(geomParams.L_d[0] / (geomParams.kappa * geomParams.h)),
        int(geomParams.L_d[1] / (geomParams.kappa * geomParams.h)),
        int(geomParams.L_d[2] / (geomParams.kappa * geomParams.h)),
        data["following_part"]["part"],
        data["following_part"]["min"],
        data["following_part"]["max"],
        data["following_part"]["particle"],
        data["following_part"]["pressure"],
        data["following_part"]["rho"],
        data["following_part"]["position"],
        data["following_part"]["velocity"],

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
        data["thermo"]["sigma"],

    };
    cout << "ThermoData initialized" << endl;


    SimulationData simParams = {

        data["simulation"]["dimension"],
        data["simulation"]["nstepT"],
        data["simulation"]["nsave"],
        data["simulation"]["dt"],
        data["simulation"]["theta"],
        data["simulation"]["alpha"],
        data["simulation"]["beta"],
        data["simulation"]["alpha_st"],
        data["simulation"]["beta_adh"],
        schemeIntegration,
        data["thermo"]["u_init"],
        state_equation,
        state_initial_condition,
        data["forces"]["gravity"],
        data["forces"]["surface_tension"],
        data["forces"]["adhesion"],
        data["condition"]["print_debug"],
        0,
        0,
        0,
        
    };
    cout << "SimulationData initialized" << endl;

    /*------------------- INITIALIZATION OF VARIABLES USED -------------------*/    

    // Initialize particles
    int MP_count = 0, FP_count = 0, GP_count = 0;
    vector<double> pos;
    vector<double> type;
    

    meshcube(geomParams, simParams, pos, type, MP_count, FP_count); 
    meshPostProcess(geomParams, simParams, pos, type, GP_count);
    
    if (!checkParticleGeneration(pos, simParams))
        exit(1);
    
    int nb_tot_part = pos.size()/3;
    vector<double> mass(nb_tot_part, 0), u(3 * nb_tot_part, 0),
                   drhodt(nb_tot_part, 0), rho(nb_tot_part, 0),
                   dudt(3 * nb_tot_part, 0), p(nb_tot_part, 0),
                   c(nb_tot_part, 0), grad_sum(nb_tot_part, 0),
                   normal(3*nb_tot_part, 0.0),
                   track_particle(simParams.nb_moving_part, 0),
                   nb_neighbours(nb_tot_part, 0);

    vector<vector<int>> cell_matrix(geomParams.Nx * geomParams.Ny * geomParams.Nz);
    vector<vector<double>> gradW(nb_tot_part),
                           W(nb_tot_part),
                           gradW_st(nb_tot_part),
                           W_st(nb_tot_part),
                           viscosity(nb_tot_part);

    int nb_sector = (simParams.dimension == 3)? 32: 8;
    vector<int> neighbours(100*nb_tot_part), 
                free_surface(nb_sector*nb_tot_part, 0);
  
    // Variables defined to used "export.cpp"
    map<string, vector<double> *> scalars;
    map<string, vector<double> *> vectors;
    scalars["type"] = &type;
    //scalars["track_particle"] = &track_particle;
    scalars["nb_neighbours"] = &nb_neighbours;
    scalars["mass"] = &mass;
    scalars["rho"] = &rho;
    scalars["p"] = &p;
    scalars["drhodt"] = &drhodt;
    scalars["grad_sum"] = &grad_sum;
    vectors["position"] = &pos;
    vectors["u"] = &u;
    vectors["dudt"] = &dudt;
    vectors["normal"] = &normal;

    printParams(geomParams, thermoParams, simParams,
                state_equation, state_initial_condition,
                MP_count, FP_count, GP_count, nb_tot_part);

    /*---------------------------- SPH ALGORITHM  ----------------------------*/

    // Initialize variables of the problem
    initRho(thermoParams, simParams, pos, rho);
    //printArray(rho,rho.size(),"rho");
    initMass(geomParams, simParams, rho, mass);
    //printArray(mass, mass.size(),"mass");
    initVelocity(thermoParams, simParams, u);
    //printArray(u, u.size(),"u");
    setPressure(geomParams, thermoParams, simParams, p, rho); 
    
    setSpeedOfSound(geomParams, thermoParams, simParams, c, rho);
    initKernelCoef(geomParams, simParams);

    
    auto t_mid = chrono::high_resolution_clock::now();
    double sim_time = 0.0;
    for (int t = 0; t < simParams.nstepT; t++){

        simParams.t = t;

        // Apply the linked-list algorithm
        sortedList(geomParams, simParams, cell_matrix, neighbours,
                   gradW, W, viscosity, nb_neighbours, type, pos, free_surface);
    
        // Compute ∇_a(W_ab) for all particles
        computeGradW(geomParams, simParams, gradW, W, neighbours, nb_neighbours, pos);
        
        // Update density, velocity and position (Euler explicit or RK22 scheme)
        updateVariables(geomParams, thermoParams, simParams, pos, u, rho, drhodt, c, p, dudt, mass, 
                        viscosity, gradW, W, neighbours, nb_neighbours, type, track_particle);


        sim_time += simParams.dt;
        // Save data each "nsave" iterations
        if(t % simParams.nsave == 0){
            if(geomParams.post_process_do || geomParams.following_part_bool ||
               geomParams.following_part_max || geomParams.following_part_min){

                    if (geomParams.post_process_do)
                        extractData(geomParams, simParams, thermoParams, pos, p, mass, neighbours, nb_neighbours, rho);
                    if (geomParams.following_part_bool || geomParams.following_part_max || geomParams.following_part_min)
                        follow_part_data(geomParams,simParams, p, rho, pos, u);
                    writing_time(sim_time);
               }
                
                
            export_particles("../../output/sph", t, pos, scalars, vectors, false);
            writing_time(sim_time);
            
            auto t_act = chrono::high_resolution_clock::now();
            double elapsed_time = double(chrono::duration_cast<chrono::duration<double>>(t_act - t_mid).count());
            progressBar(double(t)/double(simParams.nstepT), elapsed_time);    
        }
        

        // Clear matrices and reset arrays to 0
        clearAllVectors(simParams, viscosity, neighbours,
                        cell_matrix, gradW, W, drhodt, dudt, track_particle);

    }


    auto t1 = chrono::high_resolution_clock::now();
    auto delta_t = chrono::duration_cast<chrono::duration<double>>(t1 - t0).count();
    cout << "\nduration: " << delta_t << "s (" << int((MP_count+FP_count)/delta_t)*simParams.nstepT<<" [Particles*Updates/s]).\n";
    cout << "Simulation done." << endl;

    return 0;
}
