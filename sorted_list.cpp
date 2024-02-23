#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
#include<cstdlib>
#include <random>
#include <list>
#include <unordered_map>

using namespace std;


// Définition de la structure pour représenter un triplet de valeurs doubles
struct Triplet {
    double x;
    double y;
    double z;

    // Constructeur
    Triplet(double x_val, double y_val, double z_val) : x(x_val), y(y_val), z(z_val) {}
};


// Attention !! s'assurer que rand ne renvoie jamais les mêmes indices !! sinon conflit TO DO !
void setRandomParticles(const unsigned &seed, vector<vector<vector<double>>> &particule_grid, const int &nb_particles, const double &Lx, const double &Ly, const double &Lz, vector<double> &particle_x, vector<double> &particle_y, vector<double> &particle_z){

    
    // Uniform distribution 
    mt19937 gen(seed);
    uniform_real_distribution<double> dis_x(0.0, Lx); 
    uniform_real_distribution<double> dis_y(0.0, Ly); 
    uniform_real_distribution<double> dis_z(0.0, Lz); 

    for (int i = 0; i < nb_particles; ++i) {
        
        double x = dis_x(gen); 
        double y = dis_y(gen); 
        double z = dis_z(gen); 

        particle_x.push_back(x); 
        particle_y.push_back(y); 
        particle_z.push_back(z); 
    }
}


void linkedList(vector<vector<vector<double>>> &domain, vector<double> &particle_x, vector<double> &particle_y, vector<double> &particle_z, vector<unsigned> &particle_i, vector<unsigned> &particle_j, vector<unsigned> &particle_k, 
                 vector<tuple<Triplet, vector<Triplet>>> &neighbours_list, const double &Lx, const double &Ly, const double &Lz, const int &Nx, const int &Ny, const int &Nz){


    // Sort all particles in their corresponding cell
    for (double x : particle_x){
        for (double y : particle_y){
            for (double z : particle_z){
                
                particle_i.push_back(x/(Lx/Nx));
                particle_j.push_back(y/(Ly/Ny));
                particle_k.push_back(z/(Lz/Nz));
            }
        }
    }

   
    // Find neighbours for each particle
    for (unsigned x = 0; x <= particle_x.size(); x++){
        for (unsigned y = 0; y <= particle_y.size(); y++){
            for (unsigned z = 0; z <= particle_z.size(); z++){

                Triplet triplet(particle_x[x], particle_y[y], particle_z[z]); // triplet index (x,y,z) for a given particle
                vector<Triplet> associated_list_triplet; // list of neighbours for this given particle in (x,y,z)

                // Determine in which cell the particle is
                unsigned i_cell = particle_x[x]/(Lx/Nx);
                unsigned j_cell = particle_y[y]/(Ly/Ny);
                unsigned k_cell = particle_z[z]/(Lz/Nz);

                // Iterate over all 26 adjacents cells to find neighbours
                for (unsigned i = i_cell-1; i <= i_cell + 1; i++){
                    for (unsigned j = j_cell-1; j <= i_cell + 1; j++){
                        for (unsigned k = k_cell-1; k <= k_cell + 1; k++){

                            // Arrays of neighbours' index in a given cell
                            vector<int> idx_i;
                            vector<int> idx_j;
                            vector<int> idx_k;

                            // Find the index of (neighboured) particles in the cell (i,j,k)
                            for (unsigned l = 0, m = 0, n = 0 ; l < particle_i.size() && m < particle_j.size() && n < particle_k.size(); l++, m++, n++){
                                if (particle_i[l] == i){
                                    idx_i.push_back(l);
                                }
                                if (particle_i[m] == j){
                                    idx_j.push_back(m);
                                }
                                if (particle_i[n] == k){
                                    idx_k.push_back(n);
                                }
                            }

                            // Then, one knows the neighbours' indices in (i,j,k) cell
                            for (unsigned val_i = 0, val_j = 0, val_k = 0; val_i < idx_i.size(), val_j < idx_j.size(), val_k < idx_k.size(); val_i++, val_j++, val_k++){
                                associated_list_triplet.push_back(Triplet(particle_x[idx_i[val_i]], particle_y[idx_j[val_j]], particle_z[idx_k[val_k]]));
                            }

                            // For a given particle (x,y,z), add a list of triplet cordinates that gives the neighbours cordinates
                            neighbours_list.push_back((make_tuple(triplet, associated_list_triplet)));

                            idx_i.clear();
                            idx_j.clear();
                            idx_k.clear();


                        }
                    }
                }
            } 
        }
    }
} 

