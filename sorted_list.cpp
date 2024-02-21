#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
#include<cstdlib>
#include <random>




// Attention !! s'assurer que rand ne renvoie jamais les mêmes indices !! sinon conflit TO DO !
void setRandomParticles(const unsigned &seed, std::vector<std::vector<std::vector<double>>> &particule_grid, const int &nb_particles, const double &Lx, const double &Ly, const double &Lz, std::vector<double> &particle_x, std::vector<double> &particle_y, std::vector<double> &particle_z){

    
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis_x(0.0, Lx); // Distribution uniforme pour x
    std::uniform_real_distribution<double> dis_y(0.0, Ly); // Distribution uniforme pour y
    std::uniform_real_distribution<double> dis_z(0.0, Lz); // Distribution uniforme pour z

    // Génération des positions aléatoires des particules
    for (int i = 0; i < nb_particles; ++i) {
        
        double x = dis_x(gen); // Génération de x aléatoire
        double y = dis_y(gen); // Génération de y aléatoire
        double z = dis_z(gen); // Génération de z aléatoire

        particle_x.push_back(x); // Ajout de x au tableau des positions x des particules
        particle_y.push_back(y); // Ajout de y au tableau des positions y des particules
        particle_z.push_back(z); // Ajout de z au tableau des positions z des particules
    }
}


void linkedList(std::vector<std::vector<std::vector<double>>> &domain, std::vector<double> &particle_x, std::vector<double> &particle_y, std::vector<double> &particle_z, std::vector<std::tuple<double, double, double>> &neighbor_list, const double &Lx, const double &Ly, const double &Lz, const int &Nx, const int &Ny, const int &Nz){

    // Iterate over cells ONLY IF it contains at least one particle
     for (double x : particle_x){
        for (double y : particle_y){
            for (double z : particle_z){
                
                unsigned i_cell = x/(Lx/Nx);
                unsigned j_cell = y/(Ly/Ny);
                unsigned k_cell = z/(Lz/Nz);

                for (unsigned i = i_cell-1; i <= i_cell + 1; i++){
                    for (unsigned j = j_cell-1; j <= i_cell + 1; j++)}
                        for (unsigned k = k_cell-1; k <= k_cell + 1; k++){

                            /* IF VOISIN -> neighbor_list.push_back(std::make_tuple(i*(Lx/Nx), j*(Ly/Ny), k*(Lz/Nz))) */

                        }
                    }
                } 
 // 
            }
        }

