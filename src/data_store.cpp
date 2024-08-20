#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "structure.h"
#include "Kernel.h"

#include <omp.h>

using namespace std;
namespace fs = filesystem;

void delete_csvfile(string outputFile_part, SimulationData& simParams){

    if (fs::exists(outputFile_part) && simParams.t==0) {
        // Supprime le répertoire et tout son contenu
        fs::remove_all(outputFile_part);
    }
}


void extractData(GeomData &geomParams,  
                 SimulationData &simParams,
                 ThermoData &thermoParams,
                 string scheme_integration,
                 string name_file,
                 vector<double> &pos,  
                 vector<double> &p, 
                 vector<double> &mass,
                 vector<double> &u,
                 vector<int> &neighbours,
                 vector<double> &nb_neighbours,
                 vector<double> rho){

    string outputDir = "../../output/";
    if (!fs::exists(outputDir)){
        fs::create_directories(outputDir);
        cout << "Create p.csv and rho.csv and u_xyz.csv files" <<endl;
    }
    outputDir = outputDir + "/csv";
    if (!fs::exists(outputDir)){
        fs::create_directories(outputDir);
        cout << "Create following particle files" <<endl;
    }
    outputDir = outputDir +"/" + scheme_integration + "_" + name_file + "_";
    string outputFile_rho = outputDir +"rho.csv";
    string outputFile_p = outputDir + "p.csv";
    string outputFile_u_x = outputDir + "u_x.csv";
    string outputFile_u_y = outputDir + "u_y.csv";
    string outputFile_u_z = outputDir + "u_z.csv";

    // Create folder if not existing
    delete_csvfile(outputFile_rho, simParams);
    delete_csvfile(outputFile_p, simParams);
    delete_csvfile(outputFile_u_x, simParams);
    delete_csvfile(outputFile_u_y, simParams);
    delete_csvfile(outputFile_u_z, simParams);
    

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
    oss << fixed << setprecision(8) << data_part;
    

    output_name <<  data_part << "\n";

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
    oss << fixed << setprecision(8) << max;
    

    output_name << max << "\n";

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
    oss << fixed << setprecision(8) << min;
    

    output_name <<  min << "\n";

    output_name.close();
}



void follow_part_data(GeomData &geomParams,
                      SimulationData& simParams,
                      string name_file,
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

    outputDir = outputDir + "/csv";
    if (!fs::exists(outputDir)){
        fs::create_directories(outputDir);
        cout << "Create following particle files" <<endl;
    }
    string outputFile_part = outputDir + "/" + name_file + "_";


    if(geomParams.following_part_p){
        if(geomParams.following_part_bool){
            string outputFile_part_press = outputFile_part + std::to_string(particle) +"_pressure" + ".csv";
            delete_csvfile(outputFile_part_press, simParams);
            writing_in_file(outputFile_part_press, p, particle,  1, 0);
        }
        if(geomParams.following_part_max){
            string outputFile_part_press = outputFile_part + "max_pressure" + ".csv";
            delete_csvfile(outputFile_part_press, simParams);
            finding_max(outputFile_part_press,simParams, p, 1, 0);
        }

        if(geomParams.following_part_min){
            string outputFile_part_press = outputFile_part + "min_pressure" + ".csv";
            delete_csvfile(outputFile_part_press, simParams);
            finding_min(outputFile_part_press,simParams, p, 1, 0);
        }
        

    }

    if(geomParams.following_part_rho){
        if(geomParams.following_part_bool){
            string outputFile_part_rho = outputFile_part+ std::to_string(particle) +"_rho" + ".csv";
            delete_csvfile(outputFile_part_rho, simParams);
            writing_in_file(outputFile_part_rho, rho, particle,  1, 0);
        }
        if(geomParams.following_part_max){
            string outputFile_part_rho = outputFile_part + "max_rho" + ".csv";
            delete_csvfile(outputFile_part_rho, simParams);
            finding_max(outputFile_part_rho,simParams, rho, 1, 0);
        }

        if(geomParams.following_part_min){
            string outputFile_part_rho = outputFile_part + "min_rho" + ".csv";
            delete_csvfile(outputFile_part_rho, simParams);
            finding_min(outputFile_part_rho,simParams, rho, 1, 0);
        }
    }

    if(geomParams.following_part_pos[0]){
        
        if(geomParams.following_part_bool){
            string outputFile_part_pos_x = outputFile_part+ std::to_string(particle) +"_pos_x" + ".csv";
            delete_csvfile(outputFile_part_pos_x, simParams);
            writing_in_file(outputFile_part_pos_x, pos, particle,  0, 0);
        }
        if(geomParams.following_part_max){
            string outputFile_part_pos_x = outputFile_part + "max_pos_x" + ".csv";
            delete_csvfile(outputFile_part_pos_x, simParams);
            finding_max(outputFile_part_pos_x,simParams, pos, 0, 0);
        }

        if(geomParams.following_part_min){
            string outputFile_part_pos_x = outputFile_part + "min_pos_x" + ".csv";
            delete_csvfile(outputFile_part_pos_x, simParams);
            finding_min(outputFile_part_pos_x,simParams, pos, 0, 0);
        }
  
    }

    if(geomParams.following_part_pos[1]){
        
        if(geomParams.following_part_bool){
            string outputFile_part_pos_y = outputFile_part+ std::to_string(particle) +"_pos_y" + ".csv";
            delete_csvfile(outputFile_part_pos_y, simParams);
            writing_in_file(outputFile_part_pos_y, pos, particle,  0, 1);
        }
        if(geomParams.following_part_max){
            string outputFile_part_pos_y = outputFile_part + "max_pos_y" + ".csv";
            delete_csvfile(outputFile_part_pos_y, simParams);
            finding_max(outputFile_part_pos_y,simParams, pos, 0, 1);
        }

        if(geomParams.following_part_min){
            string outputFile_part_pos_y = outputFile_part + "min_pos_y" + ".csv";
            delete_csvfile(outputFile_part_pos_y, simParams);
            finding_min(outputFile_part_pos_y,simParams, pos, 0, 1);
        }
        
    }

    if(geomParams.following_part_pos[2]){
        if(geomParams.following_part_bool){
            string outputFile_part_pos_z = outputFile_part+ std::to_string(particle) +"_pos_z" + ".csv";
            delete_csvfile(outputFile_part_pos_z, simParams);
            writing_in_file(outputFile_part_pos_z, pos, particle,  0, 2);
        }
        if(geomParams.following_part_max){
            string outputFile_part_pos_z = outputFile_part + "max_pos_z" + ".csv";
            delete_csvfile(outputFile_part_pos_z, simParams);
            finding_max(outputFile_part_pos_z,simParams, pos, 0, 2);
        }

        if(geomParams.following_part_min){
            string outputFile_part_pos_z = outputFile_part + "min_pos_z" + ".csv";
            delete_csvfile(outputFile_part_pos_z, simParams);
            finding_min(outputFile_part_pos_z,simParams, pos, 0, 2);
        }
        
        
    }

    if(geomParams.following_part_u[0]){
        if(geomParams.following_part_bool){
            string outputFile_part_u_x = outputFile_part+ std::to_string(particle) +"_u_x" + ".csv";
            delete_csvfile(outputFile_part_u_x, simParams);
            writing_in_file(outputFile_part_u_x, u, particle,  0, 0);
        }
        if(geomParams.following_part_max){
            string outputFile_part_u_x = outputFile_part + "max_u_x" + ".csv";
            delete_csvfile(outputFile_part_u_x, simParams);
            finding_max(outputFile_part_u_x,simParams, u, 0, 0);
        }

        if(geomParams.following_part_min){
            string outputFile_part_u_x = outputFile_part + "min_u_x" + ".csv";
            delete_csvfile(outputFile_part_u_x, simParams);
            finding_min(outputFile_part_u_x,simParams, u, 0, 0);
        }
        
    }

    if(geomParams.following_part_u[1]){
        if(geomParams.following_part_bool){
            string outputFile_part_u_y = outputFile_part+ std::to_string(particle) +"_u_y" + ".csv";
            delete_csvfile(outputFile_part_u_y, simParams);
            writing_in_file(outputFile_part_u_y, u, particle,  0, 1);
        }
        if(geomParams.following_part_max){
            string outputFile_part_u_y = outputFile_part + "max_u_y" + ".csv";
            delete_csvfile(outputFile_part_u_y, simParams);
            finding_max(outputFile_part_u_y,simParams, u, 0, 1);
        }

        if(geomParams.following_part_min){
            string outputFile_part_u_y = outputFile_part + "min_u_y" + ".csv";
            delete_csvfile(outputFile_part_u_y, simParams);
            finding_min(outputFile_part_u_y,simParams, u, 0, 1);
        }
        
    }

    if(geomParams.following_part_u[2]){
        if(geomParams.following_part_bool){
            string outputFile_part_u_z = outputFile_part+ std::to_string(particle) +"_u_z" + ".csv";
            delete_csvfile(outputFile_part_u_z, simParams);
            writing_in_file(outputFile_part_u_z, u, particle,  0, 2);
        }
        if(geomParams.following_part_max){
            string outputFile_part_u_z = outputFile_part + "max_u_z" + ".csv";
            delete_csvfile(outputFile_part_u_z, simParams);
            finding_max(outputFile_part_u_z,simParams, u, 0, 2);
        }

        if(geomParams.following_part_min){
            string outputFile_part_u_z = outputFile_part + "min_u_z" + ".csv";
            delete_csvfile(outputFile_part_u_z, simParams);
            finding_min(outputFile_part_u_z,simParams, u, 0, 2);
        }
    }

}

void writing_time(vector<double> vec_time, string name_file, SimulationData& simParams){

    string outputDir = "../../output";

    if (!fs::exists(outputDir)){
        fs::create_directories(outputDir);
        cout << "Create following particle files" <<endl;
    }
    

    outputDir = outputDir + "/csv";
    
    string outputFile_part = outputDir + "/" + name_file + "_time" + ".csv";

    if(fs::exists(outputFile_part)){
        fs::remove_all(outputFile_part);
    }

    ofstream output_name(outputFile_part, ios::app);

    if (!output_name.is_open()){
        cerr << "Error while opening CSV file." << endl;
        return;
    }

    for(int i = 0 ; i< int(vec_time.size()) ; i++){
        ostringstream oss;
        oss << fixed << setprecision(8) << vec_time[i];
        string data_part_str = oss.str();
        output_name  << data_part_str << "\n";
    }

    output_name.close();
}
