#include <stdio.h>
#include <string.h>
#include <vector>
#include "sorted_list.h"
#include <list>
#include "cube_modified.h"
#include <iostream>
#include <fstream>
#include <cassert>

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

    // read input data from json file given as argument

    std::ifstream inputf(argv[1]);
    json data = json::parse(inputf);

    // Print input data to screen
    std::cout << argv[1] << ":\n" << data.dump(4) << std::endl;

    // Location array for particles
    vector<double> particle_x; // x-direction
    vector<double> particle_y; // y-direction
    vector<double> particle_z; // z-direction 

    std::vector<double> o = data["o"];
    std::vector<double> L = data["L"];
    double s = data["s"];
    int nstepT = data["nstepT"];

    // Number of cells we want (in each direction)
    int Nx = 4;
    int Ny = 4; 
    int Nz = 4;

    int kappa = 1;
    int h = 4;

    // cell's index for particles
    vector<unsigned> particle_i; // x-direction 
    vector<unsigned> particle_j; // y-direction
    vector<unsigned> particle_k; // z-direction

    // Initialise random particles in the domain
    meshcube(&o[0], &L[0], s, particle_x, particle_y, particle_z);

    // Location matrix for neighbours
    int nb_particles = particle_x.size();
    vector<vector<int>> neighbours_matrix_1(nb_particles);
    vector<vector<int>> neighbours_matrix_2(nb_particles);

    // Apply the linked-list algorithm
    linkedListAlgo(particle_x, particle_y, particle_z, particle_i, particle_j, particle_k, neighbours_matrix_1, &L[0], Nx, Ny, Nz, h, kappa);

    // Apply the naive algorithm
    naiveAlgo(particle_x, particle_y, particle_z, neighbours_matrix_2, h, kappa);

    for (unsigned i = 0; i < neighbours_matrix_1.size(); i++){

        std::cout << "Neighbours of particle " << i << " : ";

        for (unsigned j = 0; i < neighbours_matrix_1[i].size(); j++){

            std::cout << neighbours_matrix_1[i][j] << " (for linked list) and " << neighbours_matrix_2[i][j] << " (for naive)";
        }

        std::cout << std::endl;
    }



}