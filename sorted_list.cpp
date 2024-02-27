#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
#include <cstdlib>
#include <random>
#include <list>
#include <unordered_map>
#include <iostream>
using namespace std;



void setParticles(const int &nb_particles, const double &Lx, const double &Ly, const double &Lz, vector<double> &particle_x, vector<double> &particle_y, vector<double> &particle_z){

    for (int i = 0; i < nb_particles; ++i) {

        double x = Lx/(i+1); 
        double y = Ly/(i+1); 
        double z = Lz/(i+1); 

        particle_x.push_back(x); 
        particle_y.push_back(y); 
        particle_z.push_back(z); 
    }

}

void linkedList( vector<double> &particle_x, vector<double> &particle_y, vector<double> &particle_z, vector<unsigned> &particle_i, vector<unsigned> &particle_j, vector<unsigned> &particle_k, 
                 vector<vector<int>> neighbours_matrix, const double &Lx, const double &Ly, const double &Lz, const int &Nx, const int &Ny, const int &Nz, const int &h){

    // Sort all particles in their corresponding cell
    for (unsigned ite = 0; ite < particle_x.size(); ite ++){
                particle_i.push_back(particle_x[ite]/(Lx/Nx));
                particle_j.push_back(particle_y[ite]/(Ly/Ny));
                particle_k.push_back(particle_z[ite]/(Lz/Nz));
    }

    // Find neighbours for each particle
    for (unsigned x = 0; x < particle_x.size(); x++){
        
        // Determine in which cell the particle is a changer peut-être
        unsigned i_cell = particle_x[x]/(Lx/Nx);
        unsigned j_cell = particle_y[x]/(Ly/Ny);
        unsigned k_cell = particle_z[x]/(Lz/Nz);

        unsigned i_inf = i_cell -1;
        unsigned i_supp = i_cell+1;
        if(i_cell == 0 ){
                i_inf = 0;
        }
        if(i_cell == particle_x.size()-1 ){
                i_supp = i_cell;
        }

        unsigned j_inf = j_cell -1;
        unsigned j_supp = j_cell+1;
        if(j_cell == 0 ){
                j_inf = 0;
        }
        if(j_cell == particle_x.size()-1 ){
                j_supp = j_cell;
        }

        unsigned k_inf = k_cell -1;
        unsigned k_supp = k_cell+1;
        if(k_cell == 0 ){
                k_inf = 0;
        }
        if(k_cell == particle_x.size()-1 ){
                k_supp = k_cell;
        }

        // Iterate over all 26 adjacents cells to find neighbours 
        
        for (unsigned i = i_inf; i <= i_supp; i++){
            for (unsigned j = j_inf; j <= j_supp; j++){
                for (unsigned k = k_inf; k <= k_supp; k++){
                    for (unsigned l = x+1; l < particle_i.size(); l++){ // l est pour savoir quand on itère sur la particule si la cellule est voisine ou non 
                        if (particle_i[l] == i){
                            if (particle_j[l] == j){
                                if (particle_k[l] == k){
                                    double rx, ry, rz, r2;
                                    rx = (particle_x[x] - particle_x[l] )*(particle_x[x] - particle_x[l] );
                                    ry = (particle_y[x] - particle_y[l] )*(particle_y[x] - particle_y[l] );
                                    rz = (particle_z[x] - particle_z[l] )*(particle_z[x] - particle_z[l] );
                                    r2 = rx + ry + rz;
                                    if(r2<= k*k*h*h){
                                        neighbours_matrix[x].push_back(l); 
                                        neighbours_matrix[l].push_back(x);// ATTENTION DEMANDER VOISIN PROF KH
                                    }            
                                }      // optimisation voisin !!
                            }
                        }        
                    }
                }
            }
        }        
    }
} 
