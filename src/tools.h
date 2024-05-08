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

void progresssBar(double pourcentage, int t, int nstepT);

void clearAllVectors(SimulationData &params,
                     vector<vector<double>> &pi_matrix,
                     vector<int> &neighbours,
                     vector<vector<int>> &cell_matrix,
                     vector<vector<double>> &gradW_matrix, 
                     vector<double> &drhodt,
                     vector<double> &dudt);


#endif // TOOLS_H
