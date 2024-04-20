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
    vector<vector<double>> matrixLong = geomParams.matrixLong;
    vector<int> vectorType = geomParams.vectorType;
    double s = geomParams.s;
    int nbpart = 0;
    for(int n=0; n<vectorType.size(); n++){
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
              vector<double> &type){

    vector<vector<double>> matrixLong = geomParams.matrixLong;
    vector<vector<double>> matrixOrig = geomParams.matrixOrig;
    vector<int> vectorType = geomParams.vectorType;
    double s = geomParams.s;
    
    for(int n =0 ; n<vectorType.size();n++){
        vector<double> &L = matrixLong[n];
        vector<double> &o = matrixOrig[n];
        int type_val = vectorType[n];
        double dx = 0;
        double dy = 0;
        double dz = 0;

        int ni = int(ceil(L[0] / s));
        if(ni != 1){
            dx = L[0] / ni;
            ++ni;
        }
        
        int nj = int(ceil(L[1] / s));
        if(nj != 1){
            dy = L[1] / nj;
            ++nj;
        }
        
        int nk = int(ceil(L[2] / s));
        if(nk != 1){
            dz = L[2] / nk;
            ++nk;
        }
        
        cout<< nk<< endl;
        
        // output
        cout << "meshing cube at o=(" << o[0] << "," << o[1] << "," << o[2] << ") ";
        cout << "of size L=(" << L[0] << "," << L[1] << "," << L[2] << ")\n";
        cout << "\tparticle spacing s=(" << dx << "," << dy << "," << dz << ") [target was s=" << s << "]\n";
        cout << "\t=> " << ni << "*" << nj << "*" << nk << " = " << ni * nj * nk << " particles to be generated\n";

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
                }
            }
        }
    }
}

    


void meshBoundary(GeomData &geomParams,
                  vector<double> &bound_arr, 
                  vector<double> &type_arr){

    /*vector<double> L = geomParams.L;
    vector<double> o = geomParams.o;
    vector<double> L_d = geomParams.L_d;
    vector<double> o_d = geomParams.o_d;
    double s = geomParams.s;
    double layer_max = geomParams.particle_layers;


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

        if (geomParams.walls_used("floor")){

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


        if (geomParams.walls_used("roof")){

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


        if (geomParams.walls_used("left_wall")){

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


        if (geomParams.walls_used("right_wall")){

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


        if (geomParams.walls_used("front_wall")){

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


        if (geomParams.walls_used("back_wall")){
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
    }*/
}

void meshPostProcess(GeomData &geomParams,
                     SimulationData &simParams,
                     vector<double> &pos, 
                     vector<double> &type){

    vector<double>& post_process_in = geomParams.post_process_in;
    vector<double>& post_process_out = geomParams.post_process_out;
    double s = geomParams.s;

    double dx = post_process_out[0] - post_process_in[0];
    double dy = post_process_out[1] - post_process_in[1];
    double dz = post_process_out[2] - post_process_in[2];

    double dist = sqrt(dx * dx + dy * dy + dz * dz);

    // Particules created between initial and last particule
    int nb_points = static_cast<int>(dist / s);

    double step_x = dx / nb_points;
    double step_y = dy / nb_points;
    double step_z = dz / nb_points;
    int count = 0;

    for (int i = 0; i < nb_points; i++) {
        count++;
        pos.push_back(post_process_in[0] + i * step_x);
        pos.push_back(post_process_in[1] + i * step_y);
        pos.push_back(post_process_in[2] + i * step_z);
        type.push_back(2.0);
    }

    cout << "nb post_pro particules : " << count << endl;
}
