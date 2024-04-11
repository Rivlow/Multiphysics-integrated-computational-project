#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "structure.h"

#include <omp.h>

using namespace std;
namespace fs = std::filesystem;


void extractData(SimulationData& params,
                 vector<double> &pos, 
                 vector<double> &u, 
                 vector<double> &dudt, 
                 vector<double> &rho, 
                 vector<double> &drhodt, 
                 vector<double> &c, 
                 vector<double> &p, 
                 vector<double> &mass){

    vector<string> data_store = params.data_store;
    int init = params.data_init;
    int end = params.data_end;

    string outputDir = "../../output";


    for (auto &var : data_store){

        string outputFile = outputDir + "/" + var +".csv";

        // Create folder if not existing
        if (!fs::exists(outputDir)) {
            fs::create_directories(outputDir);
        }

        ofstream output(outputFile, ios::app);

        if (!output.is_open()) {
            cerr << "Error while opening CSV file." << endl;
            return;
        }


        // Store the desired data 
        if (var == "pos") {
            for (int i = init; i <= end; ++i) {
                if (i == end){
                    output << pos[i];
                }
                else{
                    output << pos[i] << ",";
                }
            }
            output << endl;
        }
        else if (var == "u") {
            for (int i = init; i <= end; ++i) {
            if (i == end){
                    output << u[i];
                }
                else{
                    output << u[i] << ",";
                }
            }
            output << endl;        
        }
        else if (var == "dudt") {
            for (int i = init; i <= end; ++i) {
                if (i == end){
                    output << dudt[i];
                }
                else{
                    output << dudt[i] << ",";
                }
            }
            output << endl;
        }
        else if (var == "rho") {
            for (int i = init; i <= end; ++i) {
                if (i == end){
                    output << rho[i];
                }
                else{
                    output << rho[i] << ",";
                }
            }
            output << endl;
        }
        else if (var == "drhodt") {
            for (int i = init; i <= end; ++i) {
                if (i == end){
                    output << drhodt[i];
                }
                else{
                    output << drhodt[i] << ",";
                }
            }
            output << endl;
        }
        else if (var == "c") {
            for (int i = init; i <= end; ++i) {
                if (i == end){
                    output << c[i];
                }
                else{
                    output << c[i] << ",";
                }
            }
            output<< endl;
        }
        else if (var == "p") {
            for (int i = init; i <= end; ++i) {
                if (i == end){
                    output << p[i];
                }
                else{
                    output << p[i] << ",";
                }
            }
            output << endl;
        }
        else {
            cout << "Error : data name not found." << endl;
        }

        output.close();
    }
}