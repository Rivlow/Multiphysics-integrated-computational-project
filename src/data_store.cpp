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
                 vector<double> &u,
                 vector<int> &neighbours,
                 vector<double> &nb_neighbours,
                 vector<double> rho){

    string outputDir = "../../output";
    string outputFile_rho = outputDir + "/" + "rho" +".csv";
    string outputFile_p = outputDir + "/" + "p" +".csv";
    string outputFile_u_x = outputDir + "/" + "u_x" +".csv";
    string outputFile_u_y = outputDir + "/" + "u_y" +".csv";
    string outputFile_u_z = outputDir + "/" + "u_z" +".csv";



    // Create folder if not existing
    if (!fs::exists(outputDir)){
        fs::create_directories(outputDir);
        cout << "Create p.csv and rho.csv and u_xyz.csv files" <<endl;
    }

    ofstream output_rho(outputFile_rho, ios::app);
    ofstream output_p(outputFile_p, ios::app);
    ofstream output_u_x(outputFile_u_x, ios::app);
    ofstream output_u_y(outputFile_u_y, ios::app);
    ofstream output_u_z(outputFile_u_z, ios::app);


    if (!output_rho.is_open() || !output_p.is_open() || !output_u_x.is_open() || !output_u_y.is_open() || !output_u_z.is_open()){
        cerr << "Error while opening CSV file." << endl;
        return;
    }

    int init = simParams.nb_tot_part;
    int end = pos.size()/3;

    for (int i = init; i < end; i++){

        int neighbours_size = nb_neighbours[i];
        double rho_tot = 0, p_tot = 0, u_x_tot = 0, u_y_tot = 0, u_z_tot = 0;

        for (int idx = 0; idx < neighbours_size; idx++){

            int j = neighbours[100*i + idx];
            double r_ij = 0;

            for (int coord = 0; coord < 3; coord++)
                r_ij += (pos[3*i + coord] - pos[3*j + coord])*(pos[3*i + coord] - pos[3*j + coord]);
        
            r_ij = sqrt(r_ij);
            double W_ij = CubicSpline(r_ij, geomParams, simParams);
            double m_j = mass[j];
            double rho_j = rho[j];
            double p_j = p[j];
            double u_x_j = u[3*j + 0];
            double u_y_j = u[3*j + 1];
            double u_z_j = u[3*j + 2];

            rho_tot += m_j*W_ij;
            p_tot += p_j*(m_j/rho_j)*W_ij;
            u_x_tot += u_x_j*(m_j/rho_j)*W_ij;
            u_y_tot += u_y_j*(m_j/rho_j)*W_ij;
            u_z_tot += u_z_j*(m_j/rho_j)*W_ij;

        }      
        
        if (i == end-1){
            output_rho << rho_tot << "\n";
            output_p << p_tot << "\n";
            output_u_x << u_x_tot << "\n";
            output_u_y << u_y_tot << "\n";
            output_u_z << u_z_tot << "\n";


        }
        else{
            output_rho << rho_tot << ",";
            output_p << p_tot << ",";
            output_u_x << u_x_tot << ",";
            output_u_y << u_y_tot << ",";
            output_u_z << u_z_tot << ",";
        }
    }

    output_rho.close();
    output_p.close();

}


void writingTime(double sim_time){

    string outputDir = "../../output";

    if (!fs::exists(outputDir)){
        fs::create_directories(outputDir);
        cout << "Create following particle files" <<endl;
    }
    string outputFile_part = outputDir + "/" + "time" + ".csv";

    ofstream output_name(outputFile_part, ios::app);

    if (!output_name.is_open()){
        cerr << "Error while opening CSV file." << endl;
        return;
    }

    ostringstream oss;
    oss << fixed << setprecision(6) << sim_time;
    string data_part_str = oss.str();

    for (char& ch : data_part_str) {
        if (ch == '.') {
            ch = ',';
        }
    }

    output_name << "\"" << data_part_str << "\"" << "\n";

    output_name.close();
}