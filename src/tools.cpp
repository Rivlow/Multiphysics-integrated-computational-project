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
using json = nlohmann::json;
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
// Explicit instantiation for types
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

    for (int i = 0; i < size; ++i)
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
// Explicit instantiation for types
template void printArray<int>(vector<int> &, int, string);
template void printArray<float>(vector<float> &, int, string);
template void printArray<double>(vector<double> &, int, string);



void getKey(json data,
            string &state_equation,
            string &state_initial_condition,
            string &schemeIntegration,
            vector<string> &walls_chose){
    
    for (auto &it : data["stateEquation"].items())
    {
        if (it.value() == true)
        {
            state_equation = it.key();
        }
    };

    for (auto &it : data["schemeIntegration"].items())
    {
        if (it.value() == true)
        {
            schemeIntegration = it.key();
        }
    };

    for (auto &it : data["initialCondition"].items())
    {
        if (it.value() == true)
        {
            state_initial_condition = it.key();
        }
    };

    for (auto& wall : data["domain"]["walls_used"].items()) {
        if (wall.value().get<bool>()) {
            walls_chose.push_back(wall.key());
        }
    };
}


void createOutputFolder() {

    string outputPath = "../../output";

    if (fs::exists(outputPath)) {
        cout << "Folder 'output' already exists.\n";
        return; 
    }

    try {
        fs::create_directories(outputPath);
        cout << "Folder 'output' has been successfully created.\n";
    } catch (const exception& e) {
        cerr << "Error while creating folder 'output': " << e.what() << '\n';
    }
}


void clearOutputFiles(){

    std::string outputPath = "../../output";

    try {

        for (const auto& entry : fs::directory_iterator(outputPath)){
                if (entry.path().extension() == ".vtp" || entry.path().extension() == ".csv") {
                fs::remove(entry.path());
            }
        }
        cout << "All '.vtp' files have been successfully removed.\n";

    } catch (const exception& e) {
        cerr << "Error while removing '.vtp' files : " << e.what() << '\n';
    }

}


void clearAllVectors(SimulationData &simParams,
                     vector<vector<double>> &pi_matrix,
                     vector<vector<int>> &neighbours_matrix,
                     vector<vector<int>> &cell_matrix,
                     vector<vector<double>> &gradW_matrix, 
                     vector<double> &drhodt,
                     vector<double> &dudt){

    bool PRINT = simParams.PRINT;
    int nb_moving_part = simParams.nb_moving_part;
    int cell_size = cell_matrix.size();

    for (int i = 0; i < cell_size; i++){
        cell_matrix[i].clear();
    }

    for (int i = 0; i < nb_moving_part; i++){
        gradW_matrix[i].clear();
        for(int j = 0; j < neighbours_matrix[i].size(); j++){
            neighbours_matrix[i][j] = 0.0;
        }
        
        pi_matrix[i].clear();
        drhodt[i] = 0.0;
        for(int coord = 0 ; coord < 3 ; coord ++){
            dudt[3*i+coord] = 0.0;
        
        }
        dudt[i] = 0.0;
    }


    if (PRINT){
        cout << "clearAllVectors passed" << endl;
    }


}
