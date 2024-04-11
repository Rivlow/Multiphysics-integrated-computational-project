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

    string data_store = params.data_store;
    int init = params.data_init;
    int end = params.data_end;

    string outputDir = "../../output";
    string outputFile = outputDir + "/output.csv";

    // Vérifier si le répertoire de sortie existe, sinon le créer
    if (!fs::exists(outputDir)) {
        fs::create_directories(outputDir);
    }

    ofstream output(outputFile, ios::app); // Ouvrir le fichier CSV en mode ajout

    if (!output.is_open()) {
        cerr << "Erreur lors de l'ouverture du fichier CSV." << endl;
        return;
    }

    if (!output.is_open()) {
        cerr << "Error while opening CSV file." << endl;
        return;
    }

    if (data_store == "pos") {
        for (int i = init; i <= end; ++i) {
            output << "," << pos[i];
        }
        output << endl;
    }
    else if (data_store == "u") {
        for (int i = init; i <= end; ++i) {
            output << "," << u[i];
        }
        output << endl;        
    }
    else if (data_store == "dudt") {
        for (int i = init; i <= end; ++i) {
            output << "," << dudt[i];
        }
        output << endl;
    }
    else if (data_store == "rho") {
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
    else if (data_store == "drhodt") {
        for (int i = init; i <= end; ++i) {
            output << "," << drhodt[i];
        }
        output << endl;
    }
    else if (data_store == "c") {
        for (int i = init; i <= end; ++i) {
            output << "," << c[i];
        }
        output<< endl;
    }
    else if (data_store == "p") {
        for (int i = init; i <= end; ++i) {
            output << "," << p[i];
        }
        output << endl;
    }
    else {
        cout << "Erreur de nom de données." << endl;
    }

    output.close();
}