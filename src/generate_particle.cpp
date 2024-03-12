#include "generate_particle.h"
#include <iostream>
#include <cmath> // ceil
using namespace std;

/**
 * @brief build a cube of particles aligned with x,y,z axes.
 * @param o corner of the cube with the lowest (x,y,z) values
 * @param L edge lengths along x,y and z
 * @param s particle spacing
 * @param pos_arr pos_arritions of particles
 */

int evaluateNumberParticles(vector<double> &L, double &s){

    int ni = int(ceil(L[0] / s));
    double dx = L[0] / ni;
    ++ni;
    int nj = int(ceil(L[1] / s));
    double dy = L[1] / nj;
    ++nj;
    int nk = int(ceil(L[2] / s));
    double dz = L[2] / nk;
    ++nk;

    std::cout << "\t=> " << ni << "*" << nj << "*" << nk << " = " << ni * nj * nk << " particles to be generated\n";

    return ni*nj*nk;
}

void meshcube(vector<double> &o, vector<double> &L, double &s, std::vector<double> &pos_arr)
{
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
    //std::cout << "\t=> " << ni << "*" << nj << "*" << nk << " = " << ni * nj * nk << " particles to be generated\n";

    // memory allocation
    pos_arr.reserve(pos_arr.size() + ni * nj * nk * 3);

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
                pos_arr.push_back(x);
                pos_arr.push_back(y);
                pos_arr.push_back(z);
            }
        }
    }
}

void clearAllVectors(vector<vector<double>> & artificial_visc_matrix, vector<vector<unsigned>> &neighbours_matrix, 
                     vector<vector<unsigned>> &cell_matrix, vector<vector<double>> &gradW_matrix){

    for(size_t i = 0 ; i<artificial_visc_matrix.size(); i++ ){
                
                artificial_visc_matrix[i].clear();
    }

    for(size_t i = 0 ; i<neighbours_matrix.size(); i++ ){
        
        neighbours_matrix[i].clear();
    }

    //cout << "after clear, neighbours_matrix : " << endl;

    for(size_t i = 0 ; i<cell_matrix.size(); i++ ){
        
        cell_matrix[i].clear();
    }

    //cout << "after clear, cell_matrix : " << endl;

    for(size_t i = 0 ; i<gradW_matrix.size(); i++ ){
        
        gradW_matrix[i].clear();
    }

    //cout << "after clear, gradW_matrix : " << endl;           


}
