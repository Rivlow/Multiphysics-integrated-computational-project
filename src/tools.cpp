#include "tools.h"

#include <stdio.h>
#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <cassert>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;


template <typename T>
void printMatrix(vector<vector<T>> &matrix, size_t size, string name)
{

    cout << "------------------------------------" << endl;
    cout << "           " << name << endl;
    cout << "------------------------------------"
         << "\n"
         << endl;

    for (size_t i = 0; i <= size; ++i)
    {

        cout << "For lign " << i << " : (";
        for (size_t j = 0; j < matrix[i].size(); ++j)
        {
            cout << matrix[i][j];
            if (j != matrix[i].size() - 1)
            {
                cout << ", ";
            }
        }
        cout << ")"
             << "\n"
             << endl;
    }
}
// Explicit instantiation for types we might use
template void printMatrix<int>(vector<vector<int>> &, size_t, string);
template void printMatrix<float>(vector<vector<float>> &, size_t, string);
template void printMatrix<double>(vector<vector<double>> &, size_t, string);

template <typename T>
void printArray(vector<T> &array, size_t size, string name)
{

    cout << "------------------------------------" << endl;
    cout << "           " << name << endl;
    cout << "------------------------------------"
         << "\n"
         << endl;

    for (size_t i = 0; i <= size; ++i)
    {
        std::cout << array[i];
        if (i != array.size() - 1)
        {
            std::cout << ", ";
        }
    }
    std::cout << ")"
              << "\n"
              << endl;
}
// Explicit instantiation for types we might use
template void printArray<int>(vector<int> &, size_t, string);
template void printArray<float>(vector<float> &, size_t, string);
template void printArray<double>(vector<double> &, size_t, string);

void createOutputFolder() {

    std::string outputPath = "../../output";

    if (fs::exists(outputPath)) {
        std::cout << "Folder 'output' already exists.\n";
        return; 
    }

    try {
        fs::create_directories(outputPath);
        std::cout << "Folder 'output' has been successfully created.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error while creating folder 'output': " << e.what() << '\n';
    }
}


void clearOutputFiles(){

    std::string outputPath = "../../output";

    try {

        for (const auto& entry : fs::directory_iterator(outputPath)){
            if (entry.is_regular_file() && entry.path().extension() == ".vtp") {
                // Supprimer les fichiers avec l'extension .vtp
                fs::remove(entry.path());
            }
        }
            std::cout << "All '.vtp' files have been successfully removed.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error while removing '.vtp' files : " << e.what() << '\n';
    }

}


void clearAllVectors(vector<vector<double>> &artificial_visc_matrix,
                     vector<vector<int>> &neighbours_matrix,
                     vector<vector<int>> &cell_matrix,
                     vector<vector<double>> &gradW_matrix)
{

    for (size_t i = 0; i < artificial_visc_matrix.size(); i++)
    {
        artificial_visc_matrix[i].clear();
    }
    // artificial_visc_matrix.clear();

    for (size_t i = 0; i < neighbours_matrix.size(); i++)
    {
        neighbours_matrix[i].clear();
    }
    // neighbours_matrix.clear();

    // cout << "after clear, neighbours_matrix : " << endl;

    for (size_t i = 0; i < cell_matrix.size(); i++)
    {
        cell_matrix[i].clear();
    }
    // cell_matrix.clear();

    // cout << "after clear, cell_matrix : " << endl;

    for (size_t i = 0; i < gradW_matrix.size(); i++)
    {
        gradW_matrix[i].clear();
    }
    // gradW_matrix.clear();

    // cout << "after clear, gradW_matrix : " << endl;
}
