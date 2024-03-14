#include <stdio.h>
#include <vector>
#include <string.h>
#include <list>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

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