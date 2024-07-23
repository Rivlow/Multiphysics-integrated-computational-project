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
                 vector<int> &neighbours,
                 vector<double> &nb_neighbours,
                 vector<double> rho);

void writing_in_file(string name, 
                     vector<double> data, 
                     int particle, 
                     int scalar_or_vector, 
                     int xyz);

void finding_max(string name, 
                 vector<double> data,  
                 int scalar_or_vector, 
                 int xyz);

void finding_min(string name, 
                 vector<double> data,  
                 int scalar_or_vector, 
                 int xyz);

void follow_part_data(GeomData &geomParams,
                      vector<double> p,
                      vector<double> rho,
                      vector<double> pos,
                      vector<double> u);

void writing_time(double sim_time);