#include <iostream>
#include <cmath> // ceil

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

int evaluateNumberParticles(const SimulationData &params){

    int ni = int(ceil(params.L[0] / params.s));
    // double dx = L[0] / ni;
    ++ni;
    int nj = int(ceil(params.L[1] / params.s));
    // double dy = L[1] / nj;
    ++nj;
    int nk = int(ceil(params.L[2] / params.s));
    // double dz = L[2] / nk;
    ++nk;

    std::cout << "\t=> " << ni << "*" << nj << "*" << nk << " = " << ni * nj * nk << " particles to be generated\n";

    return ni * nj * nk;
}

void meshcube(const SimulationData &params,
              vector<double> &pos_arr,
              vector<double> &type_arr){
                
    // calculate nb of particles along each direction from target size "s"
    int ni = int(ceil(params.L[0] / params.s));
    double dx = params.L[0] / ni;
    ++ni;
    int nj = int(ceil(params.L[1] / params.s));
    double dy = params.L[1] / nj;
    ++nj;
    int nk = int(ceil(params.L[2] / params.s));
    double dz = params.L[2] / nk;
    ++nk;

    // output
    std::cout << "meshing cube at o=(" << params.o[0] << "," << params.o[1] << "," << params.o[2] << ") ";
    std::cout << "of size L=(" << params.L[0] << "," << params.L[1] << "," << params.L[2] << ")\n";
    std::cout << "\tparticle spacing s=(" << dx << "," << dy << "," << dz << ") [target was s=" << params.s << "]\n";
    // std::cout << "\t=> " << ni << "*" << nj << "*" << nk << " = " << ni * nj * nk << " particles to be generated\n";
    double layer_max = params.domainParams.particle_layers;
    // memory allocation
    pos_arr.reserve(pos_arr.size() + ni * nj * nk * 3);

    // particle generation
    for (int i = 0; i < ni; ++i)
    {
        double x = params.o[0] + (layer_max-1)*params.s/2 + i * dx;
        for (int j = 0; j < nj; ++j)
        {
            double y = params.o[1] +(layer_max-1)*params.s/2 + j * dy;
            for (int k = 0; k < nk; ++k)
            {
                double z = params.o[2] +(layer_max-1)*params.s/2+ k * dz;
                pos_arr.push_back(x);
                pos_arr.push_back(y);
                pos_arr.push_back(z);

                type_arr.push_back(1.0);
            }
        }
    }
}

void meshBoundary(SimulationData &params,
                  vector<double> &bound_arr, 
                  vector<double> &type_arr){

    vector<double> L = params.L;
    vector<double> o = params.o;
    vector<double> L_d = params.L_d;
    vector<double> o_d = params.o_d;
    double s = params.s;
    double layer_max = params.domainParams.particle_layers;


    int ni;
    double dx ;
    
    int nj ;
    double dy ;
    int nk ;
    double dz ;
    
    //cout << "ni, nj, nk = " << ni << nj << nk << endl;

    //bound_arr.reserve(bound_arr.size() + 6 * (ni * nj + (nk - 2) * nj + (ni - 2) * (nk - 2)));
    int isOdd;
    
    for (int actual_layer = layer_max; 0 < actual_layer; actual_layer--){
        int ni;
        double dx ;
        int nj ;
        double dy ;
        int nk ;
        double dz ;

        double Lx = L_d[0] + (actual_layer-1)*s;
        ni = int(ceil(Lx / s));
        dx = Lx / ni;
        ++ni;
        double Ly = L_d[1] + (actual_layer-1)*s;
        nj = int(ceil(Ly / s));
        dy = Ly / nj;
        ++nj;
        double Lz = L_d[2] + (actual_layer-1)*s;
        nk = int(ceil(Lz / s));
        dz = Lz / nk;
        ++nk;
        double ox = o_d[0] + (layer_max-actual_layer)*s/2;
        double oy = o_d[1] + (layer_max-actual_layer)*s/2;
        double oz = o_d[2] + (layer_max-actual_layer)*s/2;

        if (actual_layer == 0){
            isOdd = 0;
        }
        else{
            isOdd = (actual_layer%2 != 0) ? -1 : 0;
        }

        cout<< "isOdd = " << isOdd << endl;

        if (params.walls_used("floor")){

            for (int i = 0; i < ni ; ++i){ // along x
                double x = ox + i * dx;
                for (int j = 0 ; j < nj ; ++j){ // along y
 
                    double y = oy + j * dy;
                    if (i == 0 && j == 0){
                        cout << "(x,y) : " << "(" << x << "," << y << ")" <<endl;
                    }
                    bound_arr.push_back(x);
                    bound_arr.push_back(y);
                    bound_arr.push_back(oz);
                    type_arr.push_back(0.0);
                }
            }
        }


        if (params.walls_used("roof")){

            for (int i = 0; i < ni ; ++i){ // along x
                double x = ox + i * dx;
                for (int j = 0 ; j < nj ; ++j){ // along y
                    double y = oy + j * dy;
                    double z = oz + Lz;
                    bound_arr.push_back(x);
                    bound_arr.push_back(y);
                    bound_arr.push_back(z);
                    type_arr.push_back(0.0);
                }
            }
        }


        if (params.walls_used("left_wall")){

            for (int i = 0; i < ni; ++i){ // along x
                double x = ox + i *dx;
                for (int k = 1; k < nk-1; ++k){ // along z
                
                    double z = oz + k*dz ;
                    bound_arr.push_back(x);
                    bound_arr.push_back(oy);
                    bound_arr.push_back(z);
                    type_arr.push_back(0.0);

                }
            }
        }


        if (params.walls_used("right_wall")){

            for (int i = 0; i < ni; ++i){ // along x
                double x = ox + i*dx;
                for (int k = 1; k < nk-1; ++k){ // along z
                
                    double z = oz + k*dz;
                    double y = oy + Ly;
                    bound_arr.push_back(x);
                    bound_arr.push_back(y);
                    bound_arr.push_back(z);
                    type_arr.push_back(0.0);
                    
                }
            }
        }


        if (params.walls_used("front_wall")){

            for (int j = 1; j < nj-1; ++j){ // along y
                double y = oy + j * dy;
                for (int k = 1; k < nk-1; ++k){ // along z
                
                    double z = oz + k * dz;
                    bound_arr.push_back(ox);
                    bound_arr.push_back(y);
                    bound_arr.push_back(z);
                    type_arr.push_back(0.0);
                }
            }
        }


        if (params.walls_used("back_wall")){
            for (int j = 1; j < nj-1; ++j){ // along y
                double y = oy + j*dy ;
                for (int k = 1; k < nk-1; ++k){ // along z
                    double x = ox + Lx;
                    double z = oz + k * dz;
                    bound_arr.push_back(x);
                    bound_arr.push_back(y);
                    bound_arr.push_back(z);
                    type_arr.push_back(0.0);

                }
            }
        }
    }
}