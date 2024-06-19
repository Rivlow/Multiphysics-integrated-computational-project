#include <stdio.h>
#include <vector>
#include <omp.h>

#include "NavierStokes.h"
#include "find_neighbours.h"
#include "Kernel.h"
#include "tools.h"
#include "structure.h"

using namespace std;


void surfaceTension(SimulationData& simParams,
                    GeomData &geomParams,
                    ThermoData &thermoParams,
                    vector<double> nb_neighbours,
                    vector<int> neighbours,
                    vector<vector<double>> gradW_matrix,
                    vector<vector<double>> W_matrix,
                    vector<double> mass,
                    vector<double> rho,
                    vector<double> pos,
                    vector<double> &F_vol,
                    vector<double> type,
                    vector<double> normal_grad);

void InterfaceTrackingMath(SimulationData simParams,
                           GeomData geomParams,
                           ThermoData thermoParams,
                           vector<int> nb_neighbours,
                           vector<int> neighbours,
                           vector<double> gradW,
                           vector<double> mass,
                           vector<double> rho,
                           vector<int> type,
                           vector<double> pos,
                           vector<int> &track_particle);

void surfaceTensionImprove(SimulationData& simParams,
                           GeomData &geomParams,
                           ThermoData &thermoParams,
                           vector<int> &nb_neighbours,
                           vector<int> &neighbours,
                           vector<double> &gradW,
                           vector<double> &W,
                           vector<double> &mass,
                           vector<double> &rho,
                           vector<double> &pos,
                           vector<double> &F_vol,
                           vector<int> type,
                           vector<int> &track_particle);
