#include <stdio.h>
#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <cassert>
#include <filesystem>

#include "tools.h"
#include "structure.h"



using namespace std;
namespace fs = std::filesystem;


template <typename T>
void printMatrix(vector<vector<T>> &matrix, int size, string name)
{

    cout << "------------------------------------" << endl;
    cout << "           " << name << endl;
    cout << "------------------------------------"
         << "\n"
         << endl;

    for (int i = 0; i <= size; ++i)
    {

        cout << "For lign " << i << " : (";
        for (int j = 0; j < int(matrix[i].size()); ++j)
        {
            cout << matrix[i][j];
            if (j != int(matrix[i].size() - 1))
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
template void printMatrix<int>(vector<vector<int>> &, int, string);
template void printMatrix<float>(vector<vector<float>> &, int, string);
template void printMatrix<double>(vector<vector<double>> &, int, string);

template <typename T>
void printArray(vector<T> &array, int size, string name)
{

    cout << "------------------------------------" << endl;
    cout << "           " << name << endl;
    cout << "------------------------------------"
         << "\n"
         << endl;

    for (int i = 0; i <= size; ++i)
    {
        std::cout << array[i];
        if (i != int(array.size() - 1))
        {
            std::cout << ", ";
        }
    }
    std::cout << ")"
              << "\n"
              << endl;
}
// Explicit instantiation for types we might use
template void printArray<int>(vector<int> &, int, string);
template void printArray<float>(vector<float> &, int, string);
template void printArray<double>(vector<double> &, int, string);

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


void clearAllVectors(const SimulationData &params,
                     vector<vector<double>> &artificial_visc_matrix,
                     vector<vector<int>> &neighbours_matrix,
                     vector<vector<int>> &cell_matrix,
                     vector<vector<double>> &gradW_matrix, 
                     vector<double> &drhodt_array,
                     vector<double> &dudt_array){

    bool PRINT = params.PRINT;

    for (int i = 0; i < int(artificial_visc_matrix.size()); i++){
        artificial_visc_matrix[i].clear();
    }
    // cout << "after clear, artificial_visc_matrix : " << endl;

    for (int i = 0; i < int(neighbours_matrix.size()); i++){
        neighbours_matrix[i].clear();
    }
    // cout << "after clear, neighbours_matrix : " << endl;

    for (int i = 0; i < int(cell_matrix.size()); i++){
        cell_matrix[i].clear();
    }
    // cout << "after clear, cell_matrix : " << endl;

    for (int i = 0; i < int(gradW_matrix.size()); i++){
        gradW_matrix[i].clear();
    }
    // cout << "after clear, gradW_matrix : " << endl;

    for(int i = 0 ; i < int(drhodt_array.size()); i ++ ){
        drhodt_array[i] = 0.0;
    }
    // cout << "after reset, drhodt_array : " << endl;

    for(int i = 0 ; i < int(dudt_array.size()); i ++ ){
        dudt_array[i] = 0.0;
    }
    // cout << "after reset, dudt_array : " << endl;


    if (PRINT){
        cout << "clearAllVectors passed" << endl;
    }


}
