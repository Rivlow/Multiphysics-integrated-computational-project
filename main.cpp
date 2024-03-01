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

    std::vector<double> o = data["o"];
    std::vector<double> L = data["L"];
    double s = data["s"];
    int nstepT = data["nstepT"];

    // Number of cells we want (in each direction)
    unsigned Nx, Ny, Nz ;

    int kappa = data["kappa"];
    double h = 1.2*s;

    Nx = (int) L[0]/(kappa*h); //nombre de cellules dans x 
    Ny = (int) L[1]/(kappa*h);
    Nz = (int) L[2]/(kappa*h);

    printf("(Nx, Ny, Nz) = (%d, %d, %d) \n", Nx,Ny,Nz);

    // Initialise random particles in the domain
    vector<double> part_pos;        
    vector<vector<unsigned>> cell_pos(Nx*Ny*Nz);

    meshcube(&o[0], &L[0], s, part_pos);

    // Location matrix for neighbours
    unsigned nb_particles = part_pos.size()/3;
    printf("nb of particles = %d\n",nb_particles);
    vector<vector<unsigned>> neighbours_matrix_1(nb_particles);
    vector<vector<unsigned>> neighbours_matrix_2(nb_particles);

    // Apply the linked-list algorithm
    auto start_linked = std::chrono::high_resolution_clock::now();
    linkedListAlgo(part_pos, cell_pos, neighbours_matrix_1, &L[0], Nx, Ny, Nz, h, kappa);
    auto end_linked= std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_linked = end_linked - start_linked;
    std::cout << "Time in linked list: " << elapsed_linked.count() << " secondes." << std::endl;


    // Apply the naive algorithm
    auto start_naive = std::chrono::high_resolution_clock::now();
    naiveAlgo(part_pos, neighbours_matrix_2, h, kappa);
    auto end_naive= std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_naive = end_naive - start_naive;
    std::cout << "Time in naive algo: " << elapsed_naive.count() << " secondes. \n" << std::endl;

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
        std::cout << "} (naive)\n \n";
        
    }
}