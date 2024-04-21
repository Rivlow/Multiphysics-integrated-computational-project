#include <stdio.h>
#include <vector>
#include <omp.h>

#include "gradient.h"
#include "find_neighbours.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"

using namespace std;


void surfaceTension(SimulationData& simParams, GeomData &geomParams, vector<double> nb_neighbours,
                    vector<vector<int>> neighbours_matrix, vector<double> mass,
                    vector<double> rho, vector<double> pos, vector<double> &F_vol);