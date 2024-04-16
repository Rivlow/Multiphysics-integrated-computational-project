#include <stdio.h>
#include <vector>
#include <omp.h>

#include "gradient.h"
#include "find_neighbours.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"

using namespace std;
void surfaceParticle(SimulationData& params, vector<vector<int>> neighbours_matrix,vector<double> pos, vector<double> mass,
                     vector<double> rho, vector<vector<double>> gradW_matrix,vector<int> &surface_part);

void surfaceTension(SimulationData& params,
vector<vector<double>> gradW_matrix,
vector<vector<int>> neighbours_matrix,
vector<double> mass,
vector<double> rho,
vector<double> pos,
vector<double> &F_vol
);