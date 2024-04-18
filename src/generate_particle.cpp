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

int evaluateNumberParticles(GeomData &geomParams){


    vector<double> L = geomParams.L;
    double s = geomParams.s;

    int ni = int(ceil(L[0] / s));
    ++ni;
    int nj = int(ceil(L[1] / s));
    ++nj;
    int nk = int(ceil(L[2] / s));
    ++nk;

    std::cout << "\t=> " << ni << "*" << nj << "*" << nk << " = " << ni * nj * nk << " particles to be generated\n";

    return ni * nj * nk;
}

void meshcube(GeomData &geomParams,
              vector<double> &pos,
              vector<double> &type){

    vector<double> L = geomParams.L;
    vector<double> o = geomParams.o;
    vector<double> L_d = geomParams.L_d;
    vector<double> o_d = geomParams.o_d;
    double s = geomParams.s;
    double layer_max = geomParams.particle_layers;
                
    // calculate nb of particles along each direction from target size "s"
    int ni = int(ceil(L[0] / s));
    double dx = L[0] / ni;
    ++ni;
    int nj = int(ceil(L[1] / s));
    double dy = L[1] / nj;
    ++nj;
    int nk = int(ceil(L[2] / s));
    double dz = L[2] / nk;
    ++nk;

    // output
    std::cout << "meshing cube at o=(" << o[0] << "," << o[1] << "," << o[2] << ") ";
    std::cout << "of size L=(" << L[0] << "," << L[1] << "," << L[2] << ")\n";
    std::cout << "\tparticle spacing s=(" << dx << "," << dy << "," << dz << ") [target was s=" << s << "]\n";
    // std::cout << "\t=> " << ni << "*" << nj << "*" << nk << " = " << ni * nj * nk << " particles to be generated\n";

    // memory allocation
    pos.reserve(pos.size() + ni * nj * nk * 3);

    // particle generation
    for (int i = 0; i < ni; ++i)
    {
        double x = o[0] + (layer_max-1)*s/2 + i * dx;
        for (int j = 0; j < nj; ++j)
        {
            double y = o[1] +(layer_max-1)*s/2 + j * dy;
            for (int k = 0; k < nk; ++k)
            {
                double z = o[2] +(layer_max-1)*s/2+ k * dz;
                pos.push_back(x);
                pos.push_back(y);
                pos.push_back(z);

                type.push_back(1.0);
            }
        }
    }
}

void meshBoundary(GeomData &geomParams,
                  vector<double> &bound, 
                  vector<double> &type){

    vector<double> L = geomParams.L;
    vector<double> o = geomParams.o;
    vector<double> L_d = geomParams.L_d;
    vector<double> o_d = geomParams.o_d;
    double s = geomParams.s;
    double layer_max = geomParams.particle_layers;


    //cout << "ni, nj, nk = " << ni << nj << nk << endl;

    //bound.reserve(bound.size() + 6 * (ni * nj + (nk - 2) * nj + (ni - 2) * (nk - 2)));
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

        if (geomParams.walls_used("floor")){

            for (int i = 0; i < ni ; ++i){ // along x
                double x = ox + i * dx;
                for (int j = 0 ; j < nj ; ++j){ // along y
 
                    double y = oy + j * dy;
                    if (i == 0 && j == 0){
                        cout << "(x,y) : " << "(" << x << "," << y << ")" <<endl;
                    }
                    bound.push_back(x);
                    bound.push_back(y);
                    bound.push_back(oz);
                    type.push_back(0.0);
                }
            }
        }


        if (geomParams.walls_used("roof")){

            for (int i = 0; i < ni ; ++i){ // along x
                double x = ox + i * dx;
                for (int j = 0 ; j < nj ; ++j){ // along y
                    double y = oy + j * dy;
                    double z = oz + Lz;
                    bound.push_back(x);
                    bound.push_back(y);
                    bound.push_back(z);
                    type.push_back(0.0);
                }
            }
        }


        if (geomParams.walls_used("left_wall")){

            for (int i = 0; i < ni; ++i){ // along x
                double x = ox + i *dx;
                for (int k = 1; k < nk-1; ++k){ // along z
                
                    double z = oz + k*dz ;
                    bound.push_back(x);
                    bound.push_back(oy);
                    bound.push_back(z);
                    type.push_back(0.0);

                }
            }
        }


        if (geomParams.walls_used("right_wall")){

            for (int i = 0; i < ni; ++i){ // along x
                double x = ox + i*dx;
                for (int k = 1; k < nk-1; ++k){ // along z
                
                    double z = oz + k*dz;
                    double y = oy + Ly;
                    bound.push_back(x);
                    bound.push_back(y);
                    bound.push_back(z);
                    type.push_back(0.0);
                    
                }
            }
        }


        if (geomParams.walls_used("front_wall")){

            for (int j = 1; j < nj-1; ++j){ // along y
                double y = oy + j * dy;
                for (int k = 1; k < nk-1; ++k){ // along z
                
                    double z = oz + k * dz;
                    bound.push_back(ox);
                    bound.push_back(y);
                    bound.push_back(z);
                    type.push_back(0.0);
                }
            }
        }


        if (geomParams.walls_used("back_wall")){
            for (int j = 1; j < nj-1; ++j){ // along y
                double y = oy + j*dy ;
                for (int k = 1; k < nk-1; ++k){ // along z
                    double x = ox + Lx;
                    double z = oz + k * dz;
                    bound.push_back(x);
                    bound.push_back(y);
                    bound.push_back(z);
                    type.push_back(0.0);

                }
            }
        }
    }
}

void meshPostProcess(GeomData &geomParams,
                     vector<double> &pos, 
                     vector<double> &type){

    vector<double>& post_process_in = geomParams.post_process_in;
    vector<double>& post_process_out = geomParams.post_process_out;
    double s = geomParams.s;

    double dx = post_process_out[0] - post_process_in[0];
    double dy = post_process_out[1] - post_process_in[1];
    double dz = post_process_out[2] - post_process_in[2];

    double dist = sqrt(dx * dx + dy * dy + dz * dz);

    // Nb particules created between initial and last particule
    int nb_points = static_cast<int>(dist / s);

    double step_x = dx / nb_points;
    double step_y = dy / nb_points;
    double step_z = dz / nb_points;

    for (int i = 0; i < nb_points; i++) {
        pos.push_back(post_process_in[0] + i * step_x);
        pos.push_back(post_process_in[1] + i * step_y);
        pos.push_back(post_process_in[2] + i * step_z);
        type.push_back(2.0);
    }




}