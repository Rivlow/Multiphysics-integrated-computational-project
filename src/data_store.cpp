#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "structure.h"
#include "Kernel.h"

#include <omp.h>

using namespace std;
namespace fs = std::filesystem;


void extractData(GeomData &geomParams,  
                 SimulationData& simParams,
                 ThermoData& thermoParams,
                 vector<double> &pos,  
                 vector<double> &p, 
                 vector<double> &mass,
                 vector<vector<int>> &neighbours_matrix){


    string outputDir = "../../output";
    string outputFile_rho = outputDir + "/" + "rho" +".csv";
    string outputFile_p = outputDir + "/" + "p" +".csv";

    // Create folder if not existing
    if (!fs::exists(outputDir)) fs::create_directories(outputDir);

    ofstream output_rho(outputFile_rho, ios::app);
    ofstream output_p(outputFile_p, ios::app);

    if (!output_rho.is_open() || !output_p.is_open()) {
        cerr << "Error while opening CSV file." << endl;
        return;
    }

    int init = simParams.nb_part;
    int end = pos.size()/3;

    for (int n = init; n <= end; n++) {

        vector<int> &neighbours = neighbours_matrix[n];
        int neighbours_size = neighbours.size();
        double rho_tot = 0, p_tot = 0;

        for (int idx = 0; idx < neighbours_size; idx++){

            int i_neig = neighbours[idx];
            double r_ab = 0;
            vector<double> pos_val(3);

            for (int coord = 0; coord < 3; coord++){
                
                pos_val[coord] = pos[3 * n + coord] - pos[3 * i_neig + coord];
                r_ab += pos_val[coord]*pos_val[coord];
            }

            r_ab = sqrt(r_ab);
            double W_ab = f_cubic_spline(r_ab, geomParams.h);
            double m_b = mass[i_neig];

            rho_tot += m_b*W_ab;
        }

        string state_equation = simParams.state_equation;
        double rho_0 = thermoParams.rho_0;
        double c_0 = thermoParams.c_0;
        double R = thermoParams.R;
        double T = thermoParams.T;
        double M = thermoParams.M;
        double gamma = thermoParams.gamma;

        if (state_equation == "Ideal gaz law") p_tot = (rho_tot / rho_0 - 1) * (R * T) / M;
        else if (simParams.state_equation == "Quasi incompresible fluid"){
            double B = c_0 * c_0 * rho_0 / gamma;
            p_tot = B * (pow(rho_tot / rho_0, gamma) - 1);
        }
        else {
            cout << "Error : no state equation chosen" << endl;
            exit(1);
        }


        if (n == end){
            output_rho << rho_tot;
            output_p << p_tot;
        }
        else{
            output_rho << rho_tot << ",";
            output_p << p_tot << ",";
        }
    }

    output_rho.close();
    output_p.close();

}