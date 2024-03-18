#include <stdio.h>
#include <vector>
#include <string>


using namespace std;

void deletePreviousOutputFiles();

template<typename T>
void printMatrix(vector<vector<T>> &matrix);

template<typename T>
void printArray(vector<T> &array);

void clearAllVectors(vector<vector<double>> &artificial_visc_matrix, vector<vector<unsigned>> &neighbours_matrix, 
                     vector<vector<unsigned>> &cell_matrix, vector<vector<double>> &gradW_matrix);