#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
#include<cstdlib>


struct Particle {
    int x,y,z; // particule location (global)
    int m,n,p; // cell index (local)
    double dx, dy, dz;
    std::vector<Particle> neighbor_list;
    double val_kernel;
};



struct Cell {
    int m, n, p;
    double cell_nb_x, cell_nb_y, cell_nb_z;
    std::vector<Particle> contained_particle;

    bool isParticleListEmpty() const {
        return contained_particle.empty();
    }
};


// Attention !! s'assurer que rand ne renvoie jamais les mêmes indices !! sinon conflit TO DO !
void setRandomParticles(std::vector<std::vector<std::vector<Particle>>>& particule_grid, std::vector<std::vector<std::vector<Cell>>>& cell_grid, const int &cell_nb_x, const int &cell_nb_y, const int &cell_nb_z, const int &nb_particles, std::vector<Particle>& particles_array){

    // For each particle, create an object "Particle" and assign/store a random index x,y,z to it
    for (int i = 0; i < nb_particles; i++){

        Particle particle;
        particle.x = int(rand() % particule_grid.size());
        particle.y = int(rand() % particule_grid[0].size());
        particle.z = int(rand() % particule_grid[0][0].size());

        particule_grid[particle.x][particle.y][particle.z] = particle;
        cell_grid[particle.x/cell_nb_x][particle.y/cell_nb_y][particle.z/cell_nb_z].contained_particle.push_back(particle);

    }
}

void sliceDomainIntoCells(std::vector<std::vector<std::vector<Particle>>>& particle_grid, const double &k, const double &h, std::vector<Particle>& particles,int &cell_nb_x, int &cell_nb_y, int &cell_nb_z){ // attention !! k,h peut etre changer en int ou float si besoin

    // Set constant "seed" to generate same results at each call (only used during implementation -> will be removed later)
    std::srand(0);

    int step = k*h; // est ce que k*h peut etre un double ?? a check

    /*---------------------------------------------------------------------------*/
    /* Determine the number of cells in each dimension and the size of each cell */
    /*---------------------------------------------------------------------------*/

    int len_x =  particle_grid.size(); // First dimension
    int len_y =  particle_grid[0].size(); // Second dimension
    int len_z = particle_grid[0][0].size(); // Third dimension

    int remainder_x = len_x % step;
    int remainder_y = len_y % step;
    int remainder_z = len_z % step;

    // Determine the number of cells in each dimension
    int cell_nb_x = (len_x - remainder_x)/step;
    int cell_nb_y = (len_y - remainder_y)/step;
    int cell_nb_z = (len_z - remainder_z)/step;

    // Determine the size of each cell in each dimension
    double cell_size_x_val = step + remainder_x/cell_nb_x;
    double cell_size_y_val = step + remainder_y/cell_nb_y;
    double cell_size_z_val = step + remainder_z/cell_nb_z;

    cell_nb_x = cell_size_x_val;
    cell_nb_y = cell_size_y_val;
    cell_nb_z = cell_size_z_val;
}

/*
void assignParticleToCell(std::vector<std::vector<std::vector<Cell>>>& cell_grid, const int &cell_nb_x, const int &cell_nb_y, const int &cell_nb_z){

    for (auto& particle : particles){ // "auto" is a c++ type that automatically determine the type (useful for object)

        int cell_index_x = particle.x / cell_size_x_val;
        int cell_index_y = particle.y / cell_size_y_val;
        int cell_index_z = particle.z / cell_size_z_val;

        Cell cellule;
        cellule.m = cell_index_x;
        cellule.n = cell_index_y;
        cellule.p = cell_index_z;
        cellule.contained_particle.push_back(particle);

    }

}

*/

void linkedList(std::vector<std::vector<std::vector<Cell>>>& cell_grid, const int &m, const int &n, const int &p, const int &cell_nb_x, const int &cell_nb_y, const int &cell_nb_z){

    // Iterate over cells ONLY IF it contains at least one particle
     for (unsigned m = 0; m < cell_nb_x; m++){
        for (unsigned n = 0; n < cell_nb_y; n++){
            for (unsigned p = 0; p < cell_nb_z; p++){
                if (!cell_grid[m][n][p].isParticleListEmpty()){

                    // Iteration over all particles inside specific cell
                    for (unsigned i = 0; i < cell_grid[m][n][p].contained_particle.size(); ++i) {
                        
                        const Particle& particle = cell_grid[m][n][p].contained_particle[i]; // Access the i-th particle

                        /*TO CONTINUE : NOW HAVE TO ASSIGN NEIGHBORS TO "I-TH" PARTICLE */

                } 


                
            }
        }
    }
}

