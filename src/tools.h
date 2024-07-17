#ifndef TOOLS_H
#define TOOLS_H

#include <stdio.h>
#include <vector>
#include <string>
#include <filesystem>
#include "structure.h"


using namespace std;
using json = nlohmann::json;

template <typename T>
void printMatrix(vector<vector<T>> &matrix, int size, string name);

template <typename T>
void printArray(vector<T> &array, int size, string name);

void getKey(json data,
            string &state_equation,
            string &state_initial_condition,
            string &schemeIntegration);
    

void createOutputFolder();
void clearOutputFiles();

void progressBar(double ratio, double elapsed_time);

void clearAllVectors(SimulationData &simParams,
                     vector<vector<double>> &viscosity,
                     vector<int> &neighbours,
                     vector<vector<int>> &cell_matrix,
                     vector<vector<double>> &gradW, 
                     vector<vector<double>>&W,
                     vector<double> &drhodt,
                     vector<double> &dudt,
                     vector<double> &colour,
                     vector<double> &R,
                     vector<double> &N,
                     vector<double> &normal,
                     vector<double> &track_particle,
                     vector<double> &Kappa,
                     vector<double> &dot_produc);

void printParams(GeomData geomParams,    
                 ThermoData thermoParams,
                 SimulationData simParams,
                 string state_equation,
                 string state_initial_condition,
                 int MP_count,
                 int FP_count,
                 int GP_count,
                 int nb_tot_part);


double dotProduct(vector<double> a, vector<double> b);

double dist(vector<double> &pos, int n_1, int n_2);


#endif // TOOLS_H
