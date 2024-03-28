#include <stdio.h>
#include <vector>
#include <string>
#include <filesystem>

using namespace std;

template <typename T>
void printMatrix(vector<vector<T>> &matrix, size_t size, string name);

template <typename T>
void printArray(vector<T> &array, size_t size, string name);

void createOutputFolder();
void clearOutputFiles();

void clearAllVectors(vector<vector<double>> &artificial_visc_matrix,
                     vector<vector<int>> &neighbours_matrix,
                     vector<vector<int>> &cell_matrix,
                     vector<vector<double>> &gradW_matrix, 
                     vector<double> &drhodt_array,
                     vector<double> &dudt_array,
                     const bool PRINT);

struct SimulationData {

    int kappa;
    int nstepT;
    int nsave;
    double dt; // RB
    double h; // [m] 
    double s;
    vector<double> o;
    vector<double> L;
    vector<double> o_d;
    vector<double> L_d;
    vector<double> u_init;

    double alpha;
    double beta;

    double c_0;
    double rho_moving;
    double rho_fixed;
    double rho_0;
    double M;
    double T;
    double gamma;
    static const double R = 8.314; // [J/(K.mol)]
    static const double g = -9.81; // [m/s²]

    string state_equation; 
    string state_initial_condition;
    bool PRINT;
};

