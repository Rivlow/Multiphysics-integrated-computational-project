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





void respawnParticle(vector<double> &pos, GeomData &geomParams, SimulationData &simParams);