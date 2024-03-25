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

void deletePreviousOutputFiles()
{
    std::filesystem::path outputPath = std::filesystem::current_path().parent_path().parent_path();

    cout << outputPath << endl;

    // Vérifie si le répertoire "output" existe
    if (std::filesystem::exists(outputPath) && std::filesystem::is_directory(outputPath))
    {
        // Parcours de tous les fichiers dans le répertoire "output" et suppression
        for (const auto &entry : std::filesystem::directory_iterator(outputPath))
        {
            std::filesystem::remove(entry.path());
        }
        std::cout << "Tous les fichiers dans le répertoire 'output' ont été supprimés." << std::endl;
    }
    else
    {
        std::cerr << "Le répertoire 'output' n'existe pas ou n'est pas accessible." << std::endl;
    }
}

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
