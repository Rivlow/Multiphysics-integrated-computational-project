#include "tools.h"
#include <stdio.h>
#include <vector>
#include <string.h>
#include <list>
#include <iostream>
#include <fstream>
#include <cassert>


using namespace std;

/*
template<typename T>
void printMatrix(const vector<vector<T>> &matrix){

    for (int i = 0; i < matrix.size(); ++i){

        cout << "For idx " << i << " : (";
        for (int j = 0; j < matrix[i].size(); ++j) {
            cout << matrix[i][j];
            if (j != matrix[i].size() - 1) {
                cout << ", ";
            }
        }
        cout << ")" << endl;
    }
}
*/

/*
template<typename T>
void printArray(const vector<T> &array){

    std::cout << "For idx i : (";
    for (int i = 0; i < array.size(); ++i) {
        std::cout << array[i];
        if (i != array.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << ")" << std::endl;
}
*/

void clearAllVectors(vector<vector<double>> &artificial_visc_matrix, vector<vector<unsigned>> &neighbours_matrix, 
                     vector<vector<unsigned>> &cell_matrix, vector<vector<double>> &gradW_matrix){

    for(size_t i = 0 ; i<artificial_visc_matrix.size(); i++ ){
        artificial_visc_matrix[i].clear();
    }
    //artificial_visc_matrix.clear();    

    for(size_t i = 0 ; i<neighbours_matrix.size(); i++ ){
        neighbours_matrix[i].clear();
    }
    //neighbours_matrix.clear();    


    //cout << "after clear, neighbours_matrix : " << endl;

    for(size_t i = 0 ; i<cell_matrix.size(); i++ ){
        cell_matrix[i].clear();
    }
    //cell_matrix.clear(); 

    //cout << "after clear, cell_matrix : " << endl;

    for(size_t i = 0 ; i<gradW_matrix.size(); i++ ){
        gradW_matrix[i].clear();
    }
    //gradW_matrix.clear(); 

    //cout << "after clear, gradW_matrix : " << endl;           
}
