#ifndef GENERATE_PARTICLE_H
#define GENERATE_PARTICLE_H

#include <vector>
#include "tools.h"
#include "structure.h"




using namespace std;

int evaluateNumberParticles(GeomData &geomParams);

void meshCube(GeomData &geomParams,
              SimulationData &simParams,
              vector<double> &pos,
              vector<double> &type,
              int &MP_count,
              int &FP_count);

void meshSphere(GeomData &geomParams,
                SimulationData &simParams,
                vector<double> &pos,
                vector<double> &type,
                int &MP_count,
                int &FP_count);

void meshPostProcess(GeomData &geomParams,
                     SimulationData &simParams,
                     vector<double> &pos, 
                     vector<double> &type,
                     int &GP_count);

bool checkParticleGeneration(vector<double> pos);

#endif // GENERATE_PARTICLE_H