#include <iostream>
#include <cmath> // ceil
#include <vector>
#include <unordered_set>
#include <tuple>
#include <stdio.h>
#include <string.h>

#include "generate_particle.h"
#include "structure.h"

using namespace std;

/**
 * @brief build a cube of particles aligned with x,y,z axes.
 *
 * @param o corner of the cube with the lowest (x,y,z) values
 * @param L edge lengths along x,y and z
 * @param s particle spacing
 * @param pos_arr pos_arritions of particles
 */

int evaluateNumberParticles(GeomData &geomParams){
    
    vector<vector<double>> matrixLong = geomParams.matrix_long;
    vector<int> vectorType = geomParams.vector_type;
    double s = geomParams.s;
    int nbpart = 0;

    for(int n = 0; n < int(vectorType.size()); n++){
        if(vectorType[n]){
            vector<double> &L = matrixLong[n];    
            int ni = int(ceil(L[0] / s));
            if(ni != 1){
                ++ni;
            }
            
            int nj = int(ceil(L[1] / s));
            if(nj != 1){
                ++nj;
            }

            int nk = int(ceil(L[2] / s));
            if(nk != 1){
                ++nk;
            }
            
            nbpart += ni*nj*nk;
        }
    }
    
    return nbpart;
}

void meshcube(GeomData &geomParams,
              SimulationData &simParams,
              vector<double> &pos,
              vector<double> &type,
              int &MP_count,
              int &FP_count){

    vector<vector<double>> matrixLong = geomParams.matrix_long;
    vector<vector<double>> matrixOrig = geomParams.matrix_orig;
    vector<int> vectorType = geomParams.vector_type;
    double s = geomParams.s;
    
    for(int n = 0 ; n < int(vectorType.size()); n++){

        vector<double> &L = matrixLong[n];
        vector<double> &o = matrixOrig[n];
        int type_val = vectorType[n];
        double dx = 0;
        double dy = 0;
        double dz = 0;

        double ni = ceil(L[0] / s);
        if(ni != 1){
            dx = L[0] / ni;
            ++ni;
        }
        
        double nj = ceil(L[1] / s);
        if(nj != 1){
            dy = L[1] / nj;
            ++nj;
        }
        
        double nk = ceil(L[2] / s);
        if(nk != 1){
            dz = L[2] / nk;
            ++nk;
        }
                
        // output
        cout << "meshing cube at o=(" << o[0] << "," << o[1] << "," << o[2] << ") ";
        cout << "of size L=(" << L[0] << "," << L[1] << "," << L[2] << ")\n";
        cout << "\tparticle spacing s=(" << dx << "," << dy << "," << dz << ") [target was s=" << s << "]\n";
        cout << "\t=> " << ni << "*" << nj << "*" << nk << " = " << ni * nj * nk << " particles to be generated\n";
        cout <<"\n"<< endl;

        // memory allocation
        pos.reserve(pos.size() + ni * nj * nk * 3);

        // particle generation
        for (int i = 0; i < ni; ++i)
        {
            double x = o[0] + i * dx;
            for (int j = 0; j < nj; ++j)
            {
                double y = o[1] + j * dy;
                for (int k = 0; k < nk; ++k)
                {
                    double z = o[2] + k * dz;
                    
                    pos.push_back(x);
                    pos.push_back(y);
                    pos.push_back(z);
                    type.push_back(type_val);

                    if (type_val == 1)
                        MP_count++;
                    else if (type_val == 0)
                        FP_count++;
                }
            }
        }
    }
    
    simParams.nb_tot_part = pos.size()/3;
}


void meshPostProcess(GeomData &geomParams,
                     SimulationData &simParams,
                     vector<double> &pos, 
                     vector<double> &type,
                     int &GP_count){

    if(geomParams.post_process_do){
        vector<double>& post_process_in = geomParams.xyz_init;
        vector<double>& post_process_out = geomParams.xyz_end;
        double s = geomParams.s;

        double dx = post_process_out[0] - post_process_in[0];
        double dy = post_process_out[1] - post_process_in[1];
        double dz = post_process_out[2] - post_process_in[2];

        double dist = sqrt(dx * dx + dy * dy + dz * dz);

        // Particules created between initial and last particule
        double nb_points = (dist / s);
        nb_points += 1;

        double step_x = dx / nb_points;
        double step_y = dy / nb_points;
        double step_z = dz / nb_points;

        for (int i = 0; i < nb_points; i++) {
            GP_count++;
            pos.push_back(post_process_in[0] + i * step_x);
            pos.push_back(post_process_in[1] + i * step_y);
            pos.push_back(post_process_in[2] + i * step_z);
            type.push_back(2.0);
        }
    }
}


struct TupleHash {
    template <class T1, class T2, class T3>
    size_t operator()(const tuple<T1, T2, T3>& t) const {
        auto h1 = hash<T1>{}(get<0>(t));
        auto h2 = hash<T2>{}(get<1>(t));
        auto h3 = hash<T3>{}(get<2>(t));
        return h1 ^ h2 ^ h3;
    }
};

bool checkParticleGeneration(vector<double> pos, SimulationData &simParams){

    int nb_part = simParams.nb_tot_part; 
    cout << "number of particle is  "<< nb_part << endl;
    unordered_set<tuple<double, double, double>, TupleHash> uniqueTriplets;

    for (int n = 0; n < nb_part; n++) {
        double x = pos[3 * n];
        double y = pos[3 * n + 1];
        double z = pos[3 * n + 2];

        tuple<double, double, double> triplet = make_tuple(x, y, z);

        if (uniqueTriplets.find(triplet) != uniqueTriplets.end()) {
            cout << "Error : non unique position triplet at " << n << " : (" << x << ", " << y << ", " << z << ")\n";
            return false;
        }

        uniqueTriplets.insert(triplet);
    }

    cout << "All position triplets are unique.\n";
    return true;
}
