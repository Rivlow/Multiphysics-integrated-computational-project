#include <stdio.h>
#include <string.h>
#include <vector>
#include "sorted_list.h"
#include <list>
#include "generate_particle.h"
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

    std::cout << argv[1] << ":\n" << data.dump(4) << std::endl;     // Print input data to screen

    // Location array for particles
    vector<double> particle_x; // x-direction
    vector<double> particle_y; // y-direction
    vector<double> particle_z; // z-direction 

    std::vector<double> o = data["o"];
    std::vector<double> L = data["L"];
    double s = data["s"];
    int nstepT = data["nstepT"];

    // Number of cells we want (in each direction)
    int Nx, Ny, Nz ;
    

    int kappa = data["kappa"];
    double h = 1.2*s;

    Nx = (int) L[0]/(kappa*h); //nombre de cellules dans x 
    Ny = (int) L[1]/(kappa*h);
    Nz = (int) L[2]/(kappa*h);

    // cell's index for particles
    vector<unsigned> particle_i; // x-direction 
    vector<unsigned> particle_j; // y-direction
    vector<unsigned> particle_k; // z-direction

    // Initialise random particles in the domain
    meshcube(&o[0], &L[0], s, particle_x, particle_y, particle_z);

    // Location matrix for neighbours
    int nb_particles = particle_x.size();
    printf("nb of particles = %d\n",nb_particles);
    vector<vector<int>> neighbours_matrix_1(nb_particles);
    vector<vector<int>> neighbours_matrix_2(nb_particles);

    // Apply the linked-list algorithm
    auto start_linked = std::chrono::high_resolution_clock::now();
    linkedListAlgo(particle_x, particle_y, particle_z, particle_i, particle_j, particle_k, neighbours_matrix_1, &L[0], Nx, Ny, Nz, h, kappa);
    auto end_linked= std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_linked = end_linked - start_linked;
    std::cout << "Temps écoulé en linked list: " << elapsed_linked.count() << " secondes." << std::endl;


    // Apply the naive algorithm
    auto start_naive = std::chrono::high_resolution_clock::now();
    naiveAlgo(particle_x, particle_y, particle_z, neighbours_matrix_2, h, kappa);
    auto end_naive= std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_naive = end_naive - start_naive;
    std::cout << "Temps écoulé en naive : " << elapsed_naive.count() << " secondes." << std::endl;

    for (unsigned i = 0; i < neighbours_matrix_1.size(); i++) {
        std::cout << "Particle " << i << " : ";

        std::cout << "{";
        for (unsigned j = 0; j < neighbours_matrix_1[i].size(); j++) {
            std::cout << neighbours_matrix_1[i][j];
            if (j != neighbours_matrix_1[i].size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "} (Linked-list) VS {";
        for (unsigned j = 0; j < neighbours_matrix_2[i].size(); j++) {
            std::cout << neighbours_matrix_2[i][j];
            if (j != neighbours_matrix_2[i].size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "} (naive)\n";
    }
}