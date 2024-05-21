#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
#include <cstdlib>
#include <random>
#include <list>
#include <unordered_map>
#include <iostream>
#include <chrono>
#include <algorithm>

#include "find_neighbours.h"
#include "structure.h"

#include <omp.h>

using namespace std;

void extractData(GeomData &geomParams,  
                 SimulationData& simParams,
                 ThermoData& thermoParams,
                 vector<double> &pos,  
                 vector<double> &p, 
                 vector<double> &mass,
                 vector<double> &rho,
                 vector<int> &neighbours,
                 vector<double> &nb_neighbours);