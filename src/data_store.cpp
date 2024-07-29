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
                 vector<double> &nb_neighbours,
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


void writing_in_file(string name, 
                     vector<double> data, 
                     int particle, 
                     int scalar_or_vector, 
                     int xyz){

    ofstream output_name(name, ios::app);

    if (!output_name.is_open()){
        cerr << "Error while opening CSV file." << endl;
        return;
    }

    double data_part = 0;
    if(scalar_or_vector ){
        data_part = data[particle];
    }
    else{
        data_part = data[3*particle + xyz];
    }
    
    ostringstream oss;
    oss << fixed << setprecision(6) << data_part;
    string data_part_str = oss.str();

    for (char& ch : data_part_str) {
        if (ch == '.') {
            ch = ',';
        }
    }

    output_name << "\"" << data_part_str << "\"" << "\n";

    output_name.close();

}

void finding_max(string name, 
                 SimulationData& simParams,
                 vector<double> data,  
                 int scalar_or_vector, 
                 int xyz){
    ofstream output_name(name, ios::app);

    if (!output_name.is_open()){
        cerr << "Error while opening CSV file." << endl;
        return;
    }
    double max = 0.0;
    int multiple = (scalar_or_vector == 1) ? 1 : 3;
    
    max = data[0 + xyz];
    for(int i = 1; i < simParams.nb_moving_part; i++){
        max = (max < data[multiple*i+xyz])? data[multiple*i+xyz] : max;
    }
    ostringstream oss;
    oss << fixed << setprecision(6) << max;
    string data_part_str = oss.str();

    for (char& ch : data_part_str) {
        if (ch == '.') {
            ch = ',';
        }
    }

    output_name << "\"" << data_part_str << "\"" << "\n";

    output_name.close();
}

void finding_min(string name, 
                 SimulationData& simParams,
                 vector<double> data,  
                 int scalar_or_vector, 
                 int xyz){
    ofstream output_name(name, ios::app);

    if (!output_name.is_open()){
        cerr << "Error while opening CSV file." << endl;
        return;
    }
    
    int multiple = (scalar_or_vector == 1) ? 1 : 3;
    
    double min = data[0 + xyz];
    for(int i = 1; i < simParams.nb_moving_part; i++){
        min = (min > data[multiple*i+xyz])? data[multiple*i+xyz] : min;
        
    }
    ostringstream oss;
    oss << fixed << setprecision(6) << min;
    string data_part_str = oss.str();

    for (char& ch : data_part_str) {
        if (ch == '.') {
            ch = ',';
        }
    }

    output_name << "\"" << data_part_str << "\"" << "\n";

    output_name.close();
}

void follow_part_data(GeomData &geomParams,
                      SimulationData& simParams,
                      vector<double> p,
                      vector<double> rho,
                      vector<double> pos,
                      vector<double> u){

    string outputDir = "../../output";

    int particle = geomParams.following_part_part;
    if (!fs::exists(outputDir)){
        fs::create_directories(outputDir);
        cout << "Create following particle files" <<endl;
    }
    string outputFile_part = outputDir + "/";


    if(geomParams.following_part_p){
        if(geomParams.following_part_bool){
            string outputFile_part_press = outputFile_part+ std::to_string(particle) +"_pressure" + ".csv";
            writing_in_file(outputFile_part_press, p, particle,  1, 0);
        }
        if(geomParams.following_part_max){
            string outputFile_part_press = outputFile_part + "max_pressure" + ".csv";
            finding_max(outputFile_part_press,simParams, p, 1, 0);
        }

        if(geomParams.following_part_min){
            string outputFile_part_press = outputFile_part + "min_pressure" + ".csv";
            finding_min(outputFile_part_press,simParams, p, 1, 0);
        }
        

    }

    if(geomParams.following_part_rho){
        if(geomParams.following_part_bool){
            string outputFile_part_rho = outputFile_part+ std::to_string(particle) +"_rho" + ".csv";
            writing_in_file(outputFile_part_rho, rho, particle,  1, 0);
        }
        if(geomParams.following_part_max){
            string outputFile_part_rho = outputFile_part + "max_rho" + ".csv";
            finding_max(outputFile_part_rho,simParams, rho, 1, 0);
        }

        if(geomParams.following_part_min){
            string outputFile_part_rho = outputFile_part + "min_rho" + ".csv";
            finding_min(outputFile_part_rho,simParams, rho, 1, 0);
        }
    }

    if(geomParams.following_part_pos[0]){
        
        if(geomParams.following_part_bool){
            string outputFile_part_pos_x = outputFile_part+ std::to_string(particle) +"_pos_x" + ".csv";
            writing_in_file(outputFile_part_pos_x, pos, particle,  0, 0);
        }
        if(geomParams.following_part_max){
            string outputFile_part_pos_x = outputFile_part + "max_pos_x" + ".csv";
            finding_max(outputFile_part_pos_x,simParams, pos, 0, 0);
        }

        if(geomParams.following_part_min){
            string outputFile_part_pos_x = outputFile_part + "min_pos_x" + ".csv";
            finding_min(outputFile_part_pos_x,simParams, pos, 0, 0);
        }
  
    }

    if(geomParams.following_part_pos[1]){
        
        if(geomParams.following_part_bool){
            string outputFile_part_pos_y = outputFile_part+ std::to_string(particle) +"_pos_y" + ".csv";
            writing_in_file(outputFile_part_pos_y, pos, particle,  0, 1);
        }
        if(geomParams.following_part_max){
            string outputFile_part_pos_y = outputFile_part + "max_pos_y" + ".csv";
            finding_max(outputFile_part_pos_y,simParams, pos, 0, 1);
        }

        if(geomParams.following_part_min){
            string outputFile_part_pos_y = outputFile_part + "min_pos_y" + ".csv";
            finding_min(outputFile_part_pos_y,simParams, pos, 0, 1);
        }
        
    }

    if(geomParams.following_part_pos[2]){
        if(geomParams.following_part_bool){
            string outputFile_part_pos_z = outputFile_part+ std::to_string(particle) +"_pos_z" + ".csv";
            writing_in_file(outputFile_part_pos_z, pos, particle,  0, 2);
        }
        if(geomParams.following_part_max){
            string outputFile_part_pos_z = outputFile_part + "max_pos_z" + ".csv";
            finding_max(outputFile_part_pos_z,simParams, pos, 0, 2);
        }

        if(geomParams.following_part_min){
            string outputFile_part_pos_z = outputFile_part + "min_pos_z" + ".csv";
            finding_min(outputFile_part_pos_z,simParams, pos, 0, 2);
        }
        
        
    }

    if(geomParams.following_part_u[0]){
        if(geomParams.following_part_bool){
            string outputFile_part_u_x = outputFile_part+ std::to_string(particle) +"_u_x" + ".csv";
            writing_in_file(outputFile_part_u_x, u, particle,  0, 0);
        }
        if(geomParams.following_part_max){
            string outputFile_part_u_x = outputFile_part + "max_u_x" + ".csv";
            finding_max(outputFile_part_u_x,simParams, u, 0, 0);
        }

        if(geomParams.following_part_min){
            string outputFile_part_u_x = outputFile_part + "min_u_x" + ".csv";
            finding_min(outputFile_part_u_x,simParams, u, 0, 0);
        }
        
    }

    if(geomParams.following_part_u[1]){
        if(geomParams.following_part_bool){
            string outputFile_part_u_y = outputFile_part+ std::to_string(particle) +"_u_y" + ".csv";
            writing_in_file(outputFile_part_u_y, u, particle,  0, 1);
        }
        if(geomParams.following_part_max){
            string outputFile_part_u_y = outputFile_part + "max_u_y" + ".csv";
            finding_max(outputFile_part_u_y,simParams, u, 0, 1);
        }

        if(geomParams.following_part_min){
            string outputFile_part_u_y = outputFile_part + "min_u_y" + ".csv";
            finding_min(outputFile_part_u_y,simParams, u, 0, 1);
        }
        
    }

    if(geomParams.following_part_u[2]){
        if(geomParams.following_part_bool){
            string outputFile_part_u_z = outputFile_part+ std::to_string(particle) +"_u_z" + ".csv";
            writing_in_file(outputFile_part_u_z, u, particle,  0, 2);
        }
        if(geomParams.following_part_max){
            string outputFile_part_u_z = outputFile_part + "max_u_z" + ".csv";
            finding_max(outputFile_part_u_z,simParams, u, 0, 2);
        }

        if(geomParams.following_part_min){
            string outputFile_part_u_z = outputFile_part + "min_u_z" + ".csv";
            finding_min(outputFile_part_u_z,simParams, u, 0, 2);
        }
    }

}

void writing_time(double sim_time){

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

    

    output_name << "\"" << data_part_str << "\"" << "\n";

    output_name.close();
}