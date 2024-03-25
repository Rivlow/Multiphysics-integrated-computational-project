#include <stdio.h>
#include <vector>
#include <string>
#include <filesystem>

using namespace std;

void deletePreviousOutputFiles();

template <typename T>
void printMatrix(vector<vector<T>> &matrix, size_t size, string name);

template <typename T>
void printArray(vector<T> &array, size_t size, string name);

void createOutputFolder();
void clearOutputFiles();

void clearAllVectors(vector<vector<double>> &artificial_visc_matrix,
                     vector<vector<int>> &neighbours_matrix,
                     vector<vector<int>> &cell_matrix,
                     vector<vector<double>> &gradW_matrix);