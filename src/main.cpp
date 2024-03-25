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

#include "nlohmann/json.hpp"
using json = nlohmann::json;
using namespace std;


int main(int argc, char *argv[])
{

    /*------------- SETTING COMPILATION PARAMETERS/Folders needed ---------------------------*/

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
        if (it.value() == 1)
        {
            state_equation_chosen = it.key();
        }
    }

    for (auto &it : data["initialCondition"].items())
    {
        if (it.value() == 1)
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
    //double dt = 0.005;

    // Number of cells (in each direction)
    int Nx, Ny, Nz;
    Nx = (int)L_d[0] / (kappa * h);
    Ny = (int)L_d[1] / (kappa * h);
    Nz = (int)L_d[2] / (kappa * h);

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

    // Initialization of the problem (moving particles and fixed particles)
    meshcube(o, L, pos_array, type_arr, s);
    size_t nb_moving_part = pos_array.size() / 3;
    meshBoundary(o_d, L_d, pos_array, type_arr, s);

    int nb_tot_part = pos_array.size()/3;

    std::cout << "nb_moving_part = " << nb_moving_part << std::endl;
    std::cout << "nb_tot_part = " << nb_tot_part << std::endl;
    std::cout << "s=" << s << std::endl;
    std::cout << "kappa=" << kappa << std::endl;
    std::cout << "h=" << h << std::endl;

    vector<double> mass_array(nb_tot_part),
                   u_array(3 * nb_tot_part),
                   drhodt_array(nb_tot_part),
                   rho_array(nb_tot_part),
                   dudt_array(3 * nb_tot_part, 0.0),
                   p_array(nb_tot_part);
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
    vectors["position"] = &pos_array;
    vectors["u_array"] = &u_array;
    vectors["dudt_array"] = &dudt_array;

    cout << "len(u_array) : " << pos_array.size() << endl;
    cout << "len(pos_array) : " << u_array.size() << endl;

    initializeRho(pos_array, rho_array,
                  nb_moving_part,
                  rho_init, rho_fixed, rho_0,
                  c_0, M, g, R, T, gamma,
                  state_equation_chosen, state_initial_condition);

    if (PRINT)
    {
        cout << "initializeRho passed" << endl;
    }

    initializeMass(rho_array, mass_array, s);
    if (PRINT)
    {
        cout << "initializeMass passed" << endl;
    }

    initializeVelocity(u_array, u_init, nb_moving_part);
    if (PRINT)
    {
        cout << "initializeVelocity passed" << endl;
    }

    initializeViscosity(artificial_visc_matrix);
    if (PRINT)
    {
        cout << "initializeViscosity passed" << endl;
    }


    for (int t = 0; t < nstepT; t++)
    {
        // compute dt_max = min_a h_a/F

        // double dt_max = sqrt(kappa * h / abs(g)); // CFL condition
        // std::cout << "dt_max = " << dt_max << "    dt = " << dt << std::endl;

        vector<double> drhodt_array(nb_tot_part, 0.0), dudt_array(3 * nb_tot_part, 0.0);

        // printArray(pos_array, nb_moving_part, "pos_array");

        findNeighbours(cell_matrix, neighbours_matrix, pos_array, L_d,
                       nb_moving_part, Nx, Ny, Nz, h, kappa); // Apply the linked-list algorithm


        //naiveAlgo(neighbours_matrix, pos_array, nb_moving_part, h, kappa);

        if (PRINT)
        {
            cout << "findNeighbours passed" << endl;
        }

        // Compute ∇_a(W_ab) for all particles
        gradW(gradW_matrix, 
              neighbours_matrix, 
              pos_array,
              nb_moving_part, 
              h, Nx, Ny, Nz); 

        if (PRINT)
        {
            cout << "gradW passed" << endl;
        }

        // Compute D(rho)/Dt for all particles
        continuityEquation(neighbours_matrix, gradW_matrix, pos_array, u_array, drhodt_array,
                           rho_array, mass_array,
                           nb_moving_part, h); 
        if (PRINT)
        {
            cout << "continuityEquation passed" << endl;
        }

        // Compute pressure for all particles
        setPressure(p_array, rho_array, nb_moving_part, rho_0, c_0, R, T, M,
                    gamma, state_equation_chosen); 
        if (PRINT)
        {
            cout << "setPressure passed" << endl;
        }

        // Compute artificial viscosity Π_ab for all particles
        setArtificialViscosity(artificial_visc_matrix, neighbours_matrix, pos_array,
                               rho_array, u_array, nb_moving_part, t,
                               alpha, beta, rho_0, c_0, gamma, R, T, M, h,
                               state_equation_chosen); 
        if (PRINT)
        {
            cout << "setArtificialViscosity passed" << endl;
        }

        // Compute D(u)/Dt for all particles
        momentumEquation(neighbours_matrix, 
                         gradW_matrix,
                         artificial_visc_matrix,
                         mass_array, 
                         dudt_array,
                         rho_array, 
                         p_array,
                         nb_moving_part,
                         rho_0, c_0, 
                         gamma,
                         R, T, M,
                         g, 
                         state_equation_chosen); 

        if (PRINT)
        {
            cout << "momentumEquation passed" << endl;
        }

        // printArray(rho_array, nb_moving_part, "rho_array");
        // printArray(dudt_array, nb_moving_part, "dudt_array");

        // Update density, velocity and position for each particle (Euler explicit scheme)
        for (size_t pos = 0; pos < nb_particles; pos++)
        {

            rho_array[pos] = rho_array[pos] + dt * drhodt_array[pos];

            for (size_t cord = 0; cord < 3; cord++)
            {

                pos_array[3 * pos + cord] += dt * u_array[3 * pos + cord];
                u_array[3 * pos + cord] += dt * dudt_array[3 * pos + cord];
            }
        }

        // RB: count number of neighbours
        for(int i=0; i<nb_tot_part; i++){
            nvoisins[i] = neighbours_matrix[i].size();
        }

        clearAllVectors(artificial_visc_matrix, neighbours_matrix,
                        cell_matrix, gradW_matrix);
        if (PRINT)
        {
            cout << "clearAllVectors passed" << endl;
        }


        if(t % nsave == 0)
            export_particles("sph", t, pos_array, scalars, vectors);
    }

    std::cout << "\nSimulation done." << std::endl;
    return 0;
}
