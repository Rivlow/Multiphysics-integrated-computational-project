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
                           vector<double> nb_neighbours,
                           vector<int> neighbours,
                           vector<vector<double>> gradW,
                           vector<double> mass,
                           vector<double> rho,
                           vector<double> type,
                           vector<double> pos,
                           vector<double> &track_particle);

void surfaceTensionImprove(SimulationData& simParams,
                           GeomData &geomParams,
                           ThermoData &thermoParams,
                           vector<double> &nb_neighbours,
                           vector<int> &neighbours,
                           vector<vector<double>> &gradW,
                           vector<vector<double>> &W,
                           vector<double> &mass,
                           vector<double> &rho,
                           vector<double> &pos,
                           vector<double> type,
                           vector<double> &colour,
                           vector<double> &R,
                           vector<double> &L,
                           vector<double> &N,
                           vector<double> &normal,
                           vector<double> &acc_vol,
                           vector<double> &track_particle,
                           vector<double> &Kappa,
                           vector<double> &dot_product);
