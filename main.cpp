#include <stdio.h>
#include <string.h>
#include <vector>
#include "sorted_list.cpp"
#include <list>

using namespace std;


int main() {

    /*Dimension of the domain*/

    
    // Total lenght in each direction
    double Lx = 5;
    double Ly = 5;
    double Lz = 5;

    // Number of created particles (test)
    int nb_particles = 100;

    // Number of cells we want (in each direction)
    int Nx = 4;
    int Ny = 4; 
    int Nz = 4;

    int k = 1;
    int h = 4;

    unsigned seed = 42;

    // Location array for particles
    vector<double> particle_x; // x-direction
    vector<double> particle_y; // y-direction
    vector<double> particle_z; // z-direction 

    // cell's index for particles
    vector<unsigned> particle_i; // x-direction 
    vector<unsigned> particle_j; // y-direction
    vector<unsigned> particle_k; // z-direction

    // Location matrix for neighbours
    vector<vector<int>> neighbours_matrix(nb_particles);

    // Initialise random particles in the domain
    setRandomParticles(seed, nb_particles, Lx, Ly, Lz, particle_x, particle_y, particle_z);

    // Apply the linked-list algorithm
    linkedList( particle_x, particle_y, particle_z, particle_i, particle_j, particle_k, neighbours_matrix, Lx, Ly, Lz, Nx, Ny, Nz, h);


}