#include <stdio.h>
#include <string.h>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <cassert>

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



/**
 * @brief Dummy SPH simulation testing paraview export formats
 */

int main(int argc, char *argv[]){

    /*---------------------------- SETTING COMPILATION PARAMETERS -----------------------------------*/

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

    /*---------------------------- INPUT PARAMETERS FROM JSON FILES -----------------------------------*/


    std::ifstream inputf(argv[1]);
    json data = json::parse(inputf);

    std::cout << argv[1] << ":\n" << data.dump(4) << std::endl; // Print input data to screen

    std::vector<double> o = data["o"];
    std::vector<double> L = data["L"];
    double s = data["s"];
    int nstepT = data["nstepT"];
    std::vector<double> o_d = data["o_d"];
    std::vector<double> L_d = data["L_d"];
    
    std::string state_equation_chosen, state_initial_condition;

    for (auto& it : data["stateEquation"].items()){
        if (it.value    () == 1){
            state_equation_chosen = it.key();
        }
    }

    for (auto& it : data["initialCondition"].items()){
        if (it.value    () == 1){
            state_initial_condition = it.key();
        }
    }
    
    cout << "state equation chose : " << state_equation_chosen << " \n" <<endl; 

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
    vector<double> u_init = data["u"];

    /*---------------------------- INITIALIZATION OF VARIABLES USED -----------------------------------*/

    // Debug variable
    const bool PRINT = true;

    // Constants
    double h = 1.2*s;
    double R = 8.314; // [J/(K.mol)]
    double g = 9.81; // [m/s²]
    double dt = 0.005;

    // Number of cells (in each direction) 
    unsigned Nx, Ny, Nz ;    
    Nx = (int) L_d[0]/(kappa*h);
    Ny = (int) L_d[1]/(kappa*h);
    Nz = (int) L_d[2]/(kappa*h);

    cout << " kappa * h =" << kappa*h << endl;
    printf("(Nx, Ny, Nz) = (%d, %d, %d) \n", Nx,Ny,Nz);

    // Nb of particles along each direction from target size "s"
    int nb_particles = evaluateNumberParticles(L, s);

    // Vector used (and labelled w.r.t particles location)
    vector<double> pos_arr, bound_arr;  
    vector<vector<unsigned>> cell_matrix(Nx*Ny*Nz);
 
    // Variables defined to used "export.cpp"
    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;


  /*---------------------------- SPH ALGORITHM  -----------------------------------*/

    // Initialization of the problem (moving particles and fixed particles)
    meshcube(o, L, s, pos_arr); 
    unsigned nb_moving_part = pos_arr.size()/3;
    meshBoundary(o_d, L_d, s, pos_arr);


    for(size_t i = 0; i < 3; i++){
        o_d[i] = o_d[i] - s*0.5;
        L_d[i] = L_d[i] + s ;
        //cout << " le centre et longueur de l'axe " << i << "est "<< o_d[i] << " et  " << L_d[i] << endl;
    }

    meshBoundary(o_d, L_d, s, pos_arr);

    //cout << "nb_tot_part : " << pos_arr.size() << endl; 

    unsigned nb_tot_part = pos_arr.size();

    vector<double> mass_arr(nb_tot_part), u_arr(3*nb_tot_part), drhodt_arr(nb_tot_part), rho_arr(nb_tot_part), dudt_arr(3*nb_tot_part, 0.0), p_arr(nb_tot_part);
    vector<vector<double>> artificial_visc_matrix(nb_tot_part); 
    vector<vector<unsigned>>  neighbours_matrix(nb_tot_part), neighbours_matrix_1(nb_tot_part);

    vectors["position"] = &pos_arr;
    //vectors["velocity"] = &u_arr;

    cout << "len(u_arr) : " << pos_arr.size() << endl;
    cout << "len(pos_arr) : " << u_arr.size() << endl;

    initializeRho(nb_moving_part, pos_arr, rho_arr, rho_init, rho_fixed, rho_0, c_0, M, g, R, T, gamma, state_equation_chosen, state_initial_condition);
    if(PRINT){cout << "initializeRho passed" << endl;}

    initializeMass(rho_arr, s, mass_arr);
    if(PRINT){cout << "initializeMass passed" << endl;}
    
    initializeVelocity(nb_moving_part, u_arr, u_init);
    if(PRINT){cout << "initializeVelocity passed" << endl;}
   
    initializeViscosity(artificial_visc_matrix);
    if(PRINT){cout << "initializeViscosity passed" << endl;}

    for (int t = 0; t < nstepT; t++){


    cout << "pos_arr vals : (";
        for (size_t idx = 0; idx < pos_arr.size(); idx++){
            cout << pos_arr[idx] << ", ";
        }
        cout << ") \n" << endl; 

    



        findNeighbours(nb_moving_part, pos_arr, cell_matrix, neighbours_matrix, L_d, Nx, Ny, Nz, h, kappa); // Apply the linked-list algorithm
        if(PRINT){cout << "findNeighbours passed" << endl;}

        /*
        cout << "For iteration (gradW_matrix): " << t << endl;
        for (int i = 0; i < neighbours_matrix.size(); ++i){
            cout << "For pos " << i << " : (";
            for (int j = 0; j < neighbours_matrix[i].size(); ++j) {
                cout << neighbours_matrix[i][j];
                if (j != neighbours_matrix[i].size() - 1) {
                    cout << ", ";
                }
            }
            cout << ")" << endl;
        }
        */


        vector<vector<double>> gradW_matrix;
        gradW(nb_moving_part, pos_arr, neighbours_matrix, gradW_matrix, h, Nx, Ny, Nz); // Compute ∇_a(W_ab) for all particles
        if(PRINT){cout << "gradW passed" << endl;}



        continuityEquation(nb_moving_part, pos_arr ,u_arr, neighbours_matrix, gradW_matrix, drhodt_arr, rho_arr, mass_arr, h); // Compute D(rho)/Dt for all particles
        if(PRINT){cout << "continuityEquation passed" << endl;}
        //printArray(drhodt_arr);

        setPressure(nb_moving_part, p_arr, rho_arr, rho_0, c_0, R, T, M, gamma, state_equation_chosen); // Compute pressure for all particles
        if(PRINT){cout << "setPressure passed" << endl;}

        setArtificialViscosity(t, artificial_visc_matrix, nb_moving_part, pos_arr, neighbours_matrix, rho_arr, u_arr,
                               alpha, beta, rho_0, c_0, gamma, R, T, M, h, state_equation_chosen); // Compute artificial viscosity Π_ab for all particles
        if(PRINT){cout << "setArtificialViscosity passed" << endl;}
        //printMatrix(artificial_visc_matrix);

        
        momentumEquation(nb_moving_part, neighbours_matrix, mass_arr, gradW_matrix, dudt_arr, artificial_visc_matrix, rho_arr, 
                        rho_0, c_0, p_arr, R, T, M, gamma, g, state_equation_chosen); // Compute D(u)/Dt for all particles
        if(PRINT){cout << "momentumEquation passed" << endl;}

        // Update density, velocity and position for each particle (Euler explicit scheme)
        for(size_t pos = 0; pos < nb_particles; pos++ ){

            rho_arr[pos] = rho_arr[pos] + dt*drhodt_arr[pos];

            //cout << "for particle : " << "(";
            for (size_t cord = 0; cord < 3; cord++){

                pos_arr[3*pos+cord] += dt*u_arr[3*pos+cord];
                u_arr[3*pos+cord] += dt*dudt_arr[3*pos+cord];

                // Check boundaries (temporary)
                u_arr[3*pos+cord] = (pos_arr[3*pos+cord] < 0.0) ? 0.0 : u_arr[3*pos+cord];
                u_arr[3*pos+cord] = (pos_arr[3*pos+cord] > L_d[cord]) ? 0.0 : u_arr[3*pos+cord];
                pos_arr[3*pos+cord] = (pos_arr[3*pos+cord] < 0.0) ? 0.0 : pos_arr[3*pos+cord];
                pos_arr[3*pos+cord] = (pos_arr[3*pos+cord] > L_d[cord]) ? L_d[cord] : pos_arr[3*pos+cord];
            }
        }

        clearAllVectors(artificial_visc_matrix, neighbours_matrix, cell_matrix, gradW_matrix);
        if(PRINT){cout << "clearAllVectors passed" << endl;}
        
        export_particles("test", t, pos_arr, scalars, vectors);
        /*
        cout << "pos_arr vals : (";
        for (size_t idx = 0; idx < pos_arr.size(); idx++){
            cout << pos_arr[idx] << ", ";
        }
        cout << ") \n" << endl; 

    }
    */
    }
}


