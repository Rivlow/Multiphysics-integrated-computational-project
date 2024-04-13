#include <stdio.h>
#include <vector>
#include <omp.h>

#include "gradient.h"
#include "find_neighbours.h"
#include "Kernel_functions.h"
#include "tools.h"
#include "structure.h"

using namespace std;

void surfaceTension(SimulationData& params,
vector<vector<double>> gradW_matrix,
vector<vector<int>> neighbours_matrix,
vector<double> mass,
vector<double> rho,
vector<double> type
){
    int nb_part = type.size();
    vector<double> fprime(nb_part);
    for(int n=0; n<nb_part ; n++){

        vector<double> &gradW = gradW_matrix[n];
        int size_gradW = gradW.size();
        double f = 0 ;
        for(int idx =0; idx<size_gradW/3; idx++){
            f = f+ mass[idx];

        }

    }







}