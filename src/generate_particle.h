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
              vector<int> &type,
              int &MP_count,
              int &FP_count);

void meshPostProcess(GeomData &geomParams,
                     SimulationData &simParams,
                     vector<double> &pos, 
                     vector<int> &type,
                     int &GP_count);

#endif // GENERATE_PARTICLE_H