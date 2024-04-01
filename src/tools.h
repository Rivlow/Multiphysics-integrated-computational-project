#ifndef TOOLS_H
#define TOOLS_H

#include <stdio.h>
#include <vector>
#include <string>
#include <filesystem>
#include "structure.h"


using namespace std;


template <typename T>
void printMatrix(vector<vector<T>> &matrix, int size, string name);

template <typename T>
void printArray(vector<T> &array, int size, string name);

void createOutputFolder();
void clearOutputFiles();

void clearAllVectors(const SimulationData &params,
                     vector<vector<double>> &artificial_visc_matrix,
                     vector<vector<int>> &neighbours_matrix,
                     vector<vector<int>> &cell_matrix,
                     vector<vector<double>> &gradW_matrix, 
                     vector<double> &drhodt_array,
                     vector<double> &dudt_array);


#endif // TOOLS_H
