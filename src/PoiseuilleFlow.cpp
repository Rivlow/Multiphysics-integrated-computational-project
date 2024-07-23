#include <stdio.h>
#include <vector>
#include <omp.h>

#include "NavierStokes.h"
#include "find_neighbours.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"
#include "surface_tension.h"
using namespace std;



void respawnParticle(vector<double> &pos, GeomData &geomParams, SimulationData &simParams){

    double Lx = geomParams.L_d[0],
           Ly = geomParams.L_d[1],
           Lz = geomParams.L_d[2];

    int nb_fixed_part = simParams.nb_tot_part - simParams.nb_moving_part;

    for (int i = simParams.nb_moving_part; i < simParams.nb_tot_part; i++){

        if (pos[i] >= Lx){
            for (int coord = 0; coord < 3; coord++){
                pos[3*i + coord] = 0;
            }
        }
    }






}