#include <stdio.h>
#include <string.h>
#include <vector>
#include "sorted_list.cpp"

int main() {

    int elem_x = 10;
    int elem_y = 10;
    int elem_z = 10;
    int dx = 0.1;
    int dy = 0.1;
    int dz = 0.1;
    int nb_particles = 100;
    int k = 1;
    int h = 4;

    
    std::vector<std::vector<std::vector<Particle>>> particle_grid(elem_x, std::vector<std::vector<Particle>>(elem_y, std::vector<Particle>(elem_z))); // Grid initialisation
    std::vector<Particle> particles_array; // Array of particules 
    std::vector<Cell> cells_array; // Array of cells 

    // Initialise random particles in the domain
    setRandomParticles(particle_grid, nb_particles, particles_array);

    // Divide the whole domain into cells 
    int cell_nb_x, cell_nb_y, cell_nb_z;
    sliceDomainIntoCells(particle_grid, k, h, particles_array, cell_nb_x, cell_nb_y, cell_nb_z);

    std::vector<std::vector<std::vector<Cell>>> Cell_grid(cell_nb_x, std::vector<std::vector<Cell>>(cell_nb_y, std::vector<Cell>(cell_nb_z)));


    // Apply the linked-list algorithm
    linkedList(particles_array, cell_nb_x, cell_nb_y, cell_nb_z);


}