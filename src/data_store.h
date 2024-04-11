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

void extractData(SimulationData& params,
                 vector<double> &pos, 
                 vector<double> &u, 
                 vector<double> &dudt, 
                 vector<double> &rho, 
                 vector<double> &drhodt, 
                 vector<double> &c, 
                 vector<double> &p, 
                 vector<double> &mass);