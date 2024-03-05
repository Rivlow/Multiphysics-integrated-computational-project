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

    std::cout << argv[1] << ":\n" << data.dump(4) << std::endl;     // Print input data to screen

    std::vector<double> o = data["o"];
    std::vector<double> L = data["L"];
    double s = data["s"];
    double mass = s*s*s;
    int nstepT = data["nstepT"];
    
    std::string state_equation_chosen;

    for (auto& it : data["stateEquation"].items()){
        if (it.value    () == 1){
            state_equation_chosen = it.key();
        }
    }
    
    std::cout << "state equation chose : " << state_equation_chosen << " \n" <<endl; // pk "cout" ne marche plus ??? on a definit namespace pourtant !!!

    unsigned Nx, Ny, Nz ;    // Number of cells we want (in each direction)

    int kappa = data["kappa"];
    double h = 1.2*s;
    const double R = 8.314; // [J/(K.mol)]

    Nx = (int) L[0]/(kappa*h); //nombre de cellules dans x 
    Ny = (int) L[1]/(kappa*h);
    Nz = (int) L[2]/(kappa*h);

    printf("(Nx, Ny, Nz) = (%d, %d, %d) \n", Nx,Ny,Nz);

    vector<double> part_pos;  
    vector<vector<unsigned>> cell_pos(Nx*Ny*Nz);
    meshcube(&o[0], &L[0], s, part_pos); // Initialise random particles in the domain

    unsigned nb_particles = part_pos.size()/3;

    vector<double> drhodt_arr(nb_particles), rho_arr(nb_particles), dudt_arr(3*nb_particles), p_arr(nb_particles);
    vector<vector<unsigned>> neighbours_matrix(nb_particles); // Location matrix for neighbours
    vector<vector<double>> gradW_matrix; // Location matrix for neighbours



    /*---------------------------- SPH ALGORITHM  -----------------------------------*/

    // Apply the linked-list algorithm
    findNeighbours(part_pos, cell_pos, neighbours_matrix, &L[0], Nx, Ny, Nz, h, kappa);
    std::cout << "findNeighbours algo terminated. \n" << endl;

    gradW(part_pos, neighbours_matrix, gradW_matrix, &L[0], h, Nx, Ny, Nz);
    std::cout << "gradW algo terminated. \n" << endl;

    continuityEquation(part_pos, neighbours_matrix, gradW_matrix, drhodt_arr, rho_arr, mass, h);    
    std::cout << "continuityEquation algo terminated. \n" << endl;

    /*
    for (unsigned i = 0; i < drhodt_arr.size(); i++){
        cout << "Particle " << i << " : " << drhodt_arr[i] << " \n" << endl;
    }
    */

   momentumEquation(mass, gradW_matrix, rho_arr, p_arr, state_equation_chosen);
    


}



