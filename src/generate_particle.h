#ifndef GENERATE_PARTICLE_H
#define GENERATE_PARTICLE_H

#include <vector>
#include "tools.h"
#include "structure.h"




using namespace std;

int evaluateNumberParticles(const SimulationData &params);

void meshcube(const SimulationData &params,
              vector<double> &pos_arr,
              vector<double> &type_arr);

void meshBoundary(const SimulationData &params,
                  vector<double> &bound_arr, 
                  vector<double> &type_arr);

#endif // GENERATE_PARTICLE_H