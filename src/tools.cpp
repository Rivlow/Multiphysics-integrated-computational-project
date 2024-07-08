#include <stdio.h>
#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <cassert>
#include <filesystem>
#include <chrono>

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

    for (int i = 0; i < size; ++i)
    {

        cout << "For lign " << i << " : (";
        for (int j = 0; j < int(matrix[i].size()); ++j)
        {
            cout << setprecision(8);
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
        std::cout << setprecision(8);
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
            string &schemeIntegration){

    
    for (auto &it : data["condition"]["stateEquation"].items())
    {
        if (it.value() == true)
        {
            state_equation = it.key();
        }
    };

    for (auto &it : data["condition"]["schemeIntegration"].items())
    {
        if (it.value() == true)
        {
            schemeIntegration = it.key();
        }
    };

    for (auto &it : data["condition"]["initialCondition"].items())
    {
        if (it.value() == true)
        {
            state_initial_condition = it.key();
        }
    };

    cout <<"getKey passed" << endl;
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

void progressBar(double ratio, double elapsed_time) {

    double ratio_time = ratio * 100; 
    double remain_time = (elapsed_time/ratio_time) * (100 - ratio_time);
    const int largeurBarre = 50;
    int remplissage = ratio * largeurBarre;

    cout << "[";

    for (int i = 0; i < largeurBarre; ++i) {
        if (i < remplissage) {
            cout << "#";
        } else {
            cout << " ";
        }
    }

    cout << "] " << fixed << setprecision(2) << ratio * 100 << "% (approximated remaining time = "<< remain_time << "s )\r";
    cout.flush();
}



void clearAllVectors(SimulationData &simParams,
                     vector<double> &viscosity,
                     vector<int> &neighbours,
                     vector<vector<int>> &cell_matrix,
                     vector<double> &gradW, 
                     vector<double> &drhodt,
                     vector<double> &dudt,
                     vector<int> &track_particle){

    bool PRINT = simParams.PRINT;
    int nb_tot_part = simParams.nb_tot_part;
    int cell_size = cell_matrix.size();
    int neighbours_size = neighbours.size();


    for (int i = 0; i < cell_size; i++)
        cell_matrix[i].clear();
    

    #pragma omp parallel for
    for (int idx = 0; idx < neighbours_size; idx++){
        
        neighbours[idx] = 0;
        viscosity[idx] = 0;

        for (int coord = 0; coord < 3; coord++)
            gradW[3*(idx) + coord] = 0;
    }
    

    #pragma omp parallel for
    for (int i = 0; i < nb_tot_part; i++){

        if (i <simParams.nb_moving_part)
            track_particle[i] = 0;

        drhodt[i] = 0.0;

        for(int coord = 0 ; coord < 3 ; coord ++)
            dudt[3*i+coord] = 0.0;  
        
    }    
    
    if (PRINT) cout << "clearAllVectors passed" << endl;
}

void printParams(GeomData geomParams,    
                 ThermoData thermoParams,
                 SimulationData simParams,
                 string state_equation,
                 string state_initial_condition,
                 int MP_count,
                 int FP_count,
                 int GP_count,
                 int nb_tot_part){

    cout << "#===============================#" << endl;
    cout << "# General simulation parameters #" << endl;
    cout << "#===============================#" << "\n" << endl;
    cout << "s = " << geomParams.s << " [m]" << endl;
    cout << "kappa = " << geomParams.kappa << " [-]" << endl;
    cout << "h = " << geomParams.h << " [m]" << endl;
    cout << "nstepT = " << simParams.nstepT << " [steps]" << endl;
    cout << "nsave = " << simParams.nsave << " [steps]" << endl;
    cout << "dt = " << simParams.dt << " [s]" << endl;
    cout << "theta = " << simParams.theta << " [-]" << endl;
    cout << "alpha (artificial viscosity) = " << simParams.alpha << " [-]" << endl;
    cout << "alpha (surface tension) = " << simParams.alpha_st << " [-]" << endl;
    cout << "beta (artificial viscosity) = " << simParams.beta << " [-]" << endl;
    cout << "state equation = " << state_equation << endl;
    cout << "state initial condition = " << state_initial_condition << "\n" << endl;

    cout << "#==================#" << endl;
    cout << "# Domain variables #" << endl;
    cout << "#==================#" << "\n" << endl;
    cout << "Radius of neighbourhood (kappa * h) = " << geomParams.kappa * geomParams.h << " [m]" << endl;
    cout << "Number of cells in each direction (Nx, Ny, Nz) = (" << geomParams.Nx << ", " << geomParams.Ny << ", " << geomParams.Nz << ")" << endl;
    cout << "Number of fluid particles = " << MP_count  << " [-]" << endl;
    cout << "Number of fixed particles = " << FP_count  << " [-]" << endl;
    cout << "Number of post process particles = " << GP_count << " [-]" << endl;
    cout << "Total number of particles = " << nb_tot_part - GP_count << " [-]" << "\n" << endl;

    cout << "#==================#" << endl;
    cout << "# Thermo variables #" << endl;
    cout << "#==================#" << "\n" << endl;
    cout << "Initial speed of sound (c_0) = " << thermoParams.c_0 << " [m/s]" << endl;
    cout << "Moving particle density (rho_moving) = " << thermoParams.rho_moving << " [kg/m^3]" << endl;
    cout << "Fixed particle density (rho_fixed) = " << thermoParams.rho_fixed << " [kg/m^3]" << endl;
    cout << "Initial density (rho_0) = " << thermoParams.rho_0 << " [kg/m^3]" << endl;
    cout << "Molar mass (M) = " << thermoParams.M << " [kg/mol]" << endl;
    cout << "Heat capacity ratio (gamma) = " << thermoParams.gamma << " [-]" << endl;
    cout << "Ideal gas constant (R) = " << thermoParams.R << " [J/(mol*K)]" << endl;
    cout << "Surface tension stress (sigma) = " << thermoParams.sigma << " [N/m]" << "\n" << endl;

                
}
