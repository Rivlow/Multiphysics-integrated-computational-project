#ifndef GENERATE_PARTICLE_H
#define GENERATE_PARTICLE_H

#include <vector>
#include "tools.h"
#include "structure.h"




using namespace std;

int evaluateNumberParticles(GeomData &geomParams);

void meshcube(GeomData &geomParams,
              vector<double> &pos_arr,
              vector<double> &type_arr);

void meshBoundary(GeomData &geomParams,
                  vector<double> &bound_arr, 
                  vector<double> &type_arr);

#endif // GENERATE_PARTICLE_H