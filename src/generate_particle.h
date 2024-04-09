#ifndef GENERATE_PARTICLE_H
#define GENERATE_PARTICLE_H

#include <vector>
#include "tools.h"
#include "structure.h"




using namespace std;

int evaluateNumberParticles(SimulationData &params);

void meshcube(SimulationData &params,
              vector<double> &pos_arr,
              vector<double> &type_arr);

void meshBoundary(SimulationData &params,
                  vector<double> &bound_arr, 
                  vector<double> &type_arr);

#endif // GENERATE_PARTICLE_H