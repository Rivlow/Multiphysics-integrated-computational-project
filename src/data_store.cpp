#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "structure.h"
#include "Kernel.h"

#include <omp.h>

using namespace std;
namespace fs = filesystem;


void extractData(GeomData &geomParams,  
                 SimulationData &simParams,
                 ThermoData &thermoParams,
                 vector<double> &pos,  
                 vector<double> &p, 
                 vector<double> &mass,
                 vector<int> &neighbours,
                 vector<int> &nb_neighbours,
                 vector<double> rho){

    string outputDir = "../../output";
    string outputFile_rho = outputDir + "/" + "rho" +".csv";
    string outputFile_p = outputDir + "/" + "p" +".csv";

    // Create folder if not existing
    if (!fs::exists(outputDir)){
        fs::create_directories(outputDir);
        cout << "Create p.csv and rho.csv files" <<endl;
    }

    ofstream output_rho(outputFile_rho, ios::app);
    ofstream output_p(outputFile_p, ios::app);

    if (!output_rho.is_open() || !output_p.is_open()){
        cerr << "Error while opening CSV file." << endl;
        return;
    }

    int init = simParams.nb_tot_part;
    int end = pos.size()/3;

    for (int n = init; n < end; n++){

        int neighbours_size = nb_neighbours[n];
        double rho_tot = 0, p_tot = 0;

        for (int idx = 0; idx < neighbours_size; idx++){

            //cout << "idx : " << idx <<endl;
            int i_neig = neighbours[100*n + idx];
            double r_ab = 0;
            vector<double> pos_val(3);

            for (int coord = 0; coord < 3; coord++){
                
                pos_val[coord] = pos[3 * n + coord] - pos[3 * i_neig + coord];
                r_ab += pos_val[coord]*pos_val[coord];
            }

            r_ab = sqrt(r_ab);
            double W_ab = f_cubic_spline(r_ab, geomParams, simParams);
            double m_b = mass[i_neig];
            double rho_b = rho[i_neig];
            double p_b = p[i_neig];
            rho_tot += m_b*W_ab;
            p_tot += p_b*(m_b/rho_b)*W_ab;
        }      
        
        if (n == end-1){
            output_rho << rho_tot << "\n";
            output_p << p_tot << "\n";
        }
        else{
            output_rho << rho_tot << ",";
            output_p << p_tot << ",";
        }
    }

    output_rho.close();
    output_p.close();

}
