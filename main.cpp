#include <stdio.h>
#include <string.h>
#include <vector>
#include "sorted_list.cpp"

int main() {

    /*Dimension of the domain*/

    // element size in each direction
    double dx = 0.5; 
    double dy = 0.5;
    double dz = 0.5;
    // Number of elements in each direction
    int elem_x = 100;
    int elem_y = 100;
    int elem_z = 100;
    // Total lenght in each direction
    double Lx = dx*elem_x;
    double Ly = dx*elem_y;
    double Lz = dx*elem_z;

    // Number of created particles (test)
    int nb_particles = 100;

    // Number of cells we want (in each direction)
    int Nx = 4;
    int Ny = 4; 
    int Nz = 4;

    int k = 1;
    int h = 4;

    unsigned seed = 42;

    std::vector<std::vector<std::vector<double>>> domain(elem_x, std::vector<std::vector<double>>(elem_y, std::vector<double>(elem_z))); // Grid initialisation

    // Location array for particles
    std::vector<double> particle_x; // x-location 
    std::vector<double> particle_y; // y-location
    std::vector<double> particle_z; // z-location 

    // List of neighbors for each particle
    std::vector<std::tuple<double, double, double>> neighbor_list;


    // Initialise random particles in the domain
    setRandomParticles(seed, domain, nb_particles, Lx, Ly, Lz, particle_x, particle_y, particle_z);

    // Apply the linked-list algorithm
    linkedList(domain, particle_x, particle_y, particle_z, neighbor_list, Nx, Ny, Nz);


}