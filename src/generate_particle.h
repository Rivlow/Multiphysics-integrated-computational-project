#ifndef GENERATE_PARTICLE_H
#define GENERATE_PARTICLE_H

#include <vector>
#include "tools.h"
#include "structure.h"




using namespace std;

int evaluateNumberParticles(GeomData &geomParams);

void meshcube(GeomData &geomParams,
              SimulationData &simParams,
              vector<double> &pos,
              vector<double> &type);

void meshPostProcess(GeomData &geomParams,
                     SimulationData &simParams,
                     vector<double> &pos, 
                     vector<double> &type);

#endif // GENERATE_PARTICLE_H