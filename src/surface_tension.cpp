#include <stdio.h>
#include <vector>
#include <cmath>
#include <omp.h>

#include "gradient.h"
#include "find_neighbours.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"
#include "surface_tension.h"
#include "Kernel.h"

using namespace std;


void surfaceTension(SimulationData& simParams,
                    GeomData &geomParams,
                    ThermoData &thermoParam,
                    vector<double> nb_neighbours,
                    vector<int> neighbours,
                    vector<vector<double>> gradW_matrix,
                    vector<vector<double>> W_matrix,
                    vector<double> mass,
                    vector<double> rho,
                    vector<double> pos,
                    vector<double> &F_vol,
                    vector<double> type,
                    vector<double> normal_grad){

    vector<double> normal(3*simParams.nb_moving_part,0.0);
    #pragma omp parallel for 
    for(int n = 0; n<simParams.nb_moving_part; n++){

        vector<double> &gradW = gradW_matrix[n];
        int size_neighbours = nb_neighbours[n];
        
        for(int idx = 0; idx<size_neighbours; idx++){ 

            int i_neig = neighbours[100*n + idx]; 
            double m_j = mass[i_neig];
            double rho_j = rho[i_neig];
            
            for( int coord = 0; coord <3; coord ++){

                double grad = gradW[3*idx+coord];
                normal[3*n+coord] += geomParams.h*m_j*grad/rho_j;
            }
        }
    }
    double alpha = simParams.alpha_st;
    #pragma omp parallel for 
    for(int n = 0; n<simParams.nb_moving_part; n++){
    
        int size_neighbours = nb_neighbours[n];
        
        for(int idx = 0; idx<size_neighbours; idx++){
            
            if(type[neighbours[100*n + idx]] == 1){
                int i_neig = neighbours[100*n + idx];
                double K_ij = 2*thermoParam.rho_0/(rho[n]+rho[i_neig]);
                double r_ab = 0;
                vector<double> d_xyz(3);
                
                for (int coord = 0; coord < 3; coord++){
                    
                    d_xyz[coord] = pos[3 * n + coord] - pos[3 * i_neig + coord];
                    r_ab += d_xyz[coord]*d_xyz[coord];
                }
                
                r_ab = sqrt(r_ab);
                double W_ab = W_coh(r_ab,geomParams.kappa*geomParams.h, simParams);
                
                double m_a = mass[n];
                double m_b = mass[i_neig];
                double F_res = 0;
                
                for (int coord = 0; coord < 3; coord++){
                    
                    F_vol[3*n + coord] += -K_ij*(alpha * m_a * m_b * d_xyz[coord]*W_ab/r_ab 
                                    + alpha*(normal[3*n+coord]-normal[3*i_neig+coord]));
                
                    F_res += F_vol[3*n + coord]*F_vol[3*n + coord];
                }
                
                simParams.F_st_max = sqrt(F_res);
            }
        }     
    }  
}

/*
#include <vector>
#include <cmath>

struct Vec3 {
    double x, y, z;

    Vec3 operator+(const Vec3& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    Vec3 operator-(const Vec3& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    Vec3 operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    double dot(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    double length() const {
        return std::sqrt(x*x + y*y + z*z);
    }

    Vec3 normalize() const {
        double len = length();
        return {x / len, y / len, z / len};
    }
};

struct Particle {
    Vec3 position;
    Vec3 velocity;
    Vec3 force;
    double density;
    double pressure;
};

void detectFreeSurface(const std::vector<Particle>& particles, std::vector<bool>& isSurface, double threshold) {
    for (size_t i = 0; i < particles.size(); ++i) {
        Vec3 normal = {0, 0, 0};
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i != j) {
                Vec3 r = particles[i].position - particles[j].position;
                double dist = r.length();
                if (dist < threshold) {
                    normal = normal + (r * (1.0 / dist));
                }
            }
        }
        isSurface[i] = normal.length() > 0.5; // Exemple de seuil
    }
}

void addImaginaryParticles(const std::vector<Particle>& particles, std::vector<Particle>& imaginaryParticles, const std::vector<bool>& isSurface, double spacing) {
    for (size_t i = 0; i < particles.size(); ++i) {
        if (isSurface[i]) {
            Vec3 normal = {0, 0, 0};
            for (size_t j = 0; j < particles.size(); ++j) {
                if (i != j) {
                    Vec3 r = particles[i].position - particles[j].position;
                    double dist = r.length();
                    if (dist < spacing) {
                        normal = normal + (r * (1.0 / dist));
                    }
                }
            }
            normal = normal.normalize();
            Particle imaginaryParticle;
            imaginaryParticle.position = particles[i].position - (normal * spacing);
            imaginaryParticles.push_back(imaginaryParticle);
        }
    }
}

double computeCurvature(const Vec3& normal, double kernelRadius) {
    // Utiliser un noyau de lissage pour estimer la courbure
    // Placeholder pour le calcul de la courbure basé sur les normales
    return normal.length() / kernelRadius;
}

void computeSurfaceTension(std::vector<Particle>& particles, const std::vector<Particle>& imaginaryParticles, double sigma, double kernelRadius) {
    for (size_t i = 0; i < particles.size(); ++i) {
        Vec3 normal = {0, 0, 0};
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i != j) {
                Vec3 r = particles[i].position - particles[j].position;
                double dist = r.length();
                if (dist < kernelRadius) {
                    normal = normal + (r * (1.0 / dist));
                }
            }
        }

        double curvature = computeCurvature(normal, kernelRadius);
        Vec3 surfaceTensionForce = normal * (-sigma * curvature);
        particles[i].force = particles[i].force + surfaceTensionForce;
    }
}

int main() {
    std::vector<Particle> particles; // Initialisation des particules
    std::vector<Particle> imaginaryParticles;
    std::vector<bool> isSurface(particles.size(), false);

    double threshold = 0.1;
    double spacing = 0.02;
    double sigma = 0.0728; // Tension de surface pour l'eau à 20°C
    double kernelRadius = 0.04;

    detectFreeSurface(particles, isSurface, threshold);
    addImaginaryParticles(particles, imaginaryParticles, isSurface, spacing);
    computeSurfaceTension(particles, imaginaryParticles, sigma, kernelRadius);

    // Simulation loop
    for (int t = 0; t < 1000; ++t) {
        // Mise à jour des positions et des vitesses des particules
        // ...
    }

    return 0;
}

*/