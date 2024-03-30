

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

    // memory allocation
    pos_arr.reserve(pos_arr.size() + ni * nj * nk * 3);

    // particle generation
    for (int i = 0; i < ni; ++i)
    {
        double x = params.o[0] + i * dx;
        for (int j = 0; j < nj; ++j)
        {
            double y = params.o[1] + j * dy;
            for (int k = 0; k < nk; ++k)
            {
                double z = params.o[2] + k * dz;
                pos_arr.push_back(x);
                pos_arr.push_back(y);
                pos_arr.push_back(z);

                type_arr.push_back(1.0);
            }
        }
    }
}

void meshBoundary(const SimulationData &params,
                  vector<double> &bound_arr, 
                  vector<double> &type_arr){

    int ni = int(ceil(params.L_d[0] / params.s));
    double dx = params.L_d[0] / ni;
    ++ni;
    int nj = int(ceil(params.L_d[1] / params.s));
    double dy = params.L_d[1] / nj;
    ++nj;
    int nk = int(ceil(params.L_d[2] / params.s));
    double dz = params.L_d[2] / nk;
    ++nk;
    cout << "ni, nj, nk = " << ni << nj << nk << endl;

    bound_arr.reserve(bound_arr.size() + 6 * (ni * nj + (nk - 2) * nj + (ni - 2) * (nk - 2)));

    // Apply first layer of FP
    for (int i = 0; i < ni ; ++i) // along x
    {
        double x = params.o_d[0] + i * dx;
        for (int j = 0 ; j < nj ; ++j) // along y
        {
            double y = params.o_d[1] + j * dy;
            bound_arr.push_back(x);
            bound_arr.push_back(y);
            bound_arr.push_back(params.o_d[2]);
            type_arr.push_back(0.0);
            // bound_arr.push_back(x);
            // bound_arr.push_back(y);
            // bound_arr.push_back(L_d[2] + o_d[2]);
            // type_arr.push_back(0.0);
        }
    }

    for (int j = 0; j < nj; ++j) // along y
    {
        double y = params.o_d[1] + j * dy;
        for (int k = 1; k < nk-1; ++k) // along z
        {
         double z = params.o_d[2] + k * dz;
        bound_arr.push_back(params.o_d[0]);
        bound_arr.push_back(y);
        bound_arr.push_back(z);
        type_arr.push_back(0.0);
        bound_arr.push_back(params.L_d[0] + params.o_d[0]);
        bound_arr.push_back(y);
        bound_arr.push_back(z);
        type_arr.push_back(0.0);

        }
    }
    for (int i = 1; i < ni-1; ++i) // along x
    {
        double x = params.o_d[0] + i * dx;
        for (int k = 1; k < nk-1; ++k) // along z
        {
         double z = params.o_d[2] + k * dz;
        bound_arr.push_back(x);
        bound_arr.push_back(params.o_d[1]);
        bound_arr.push_back(z);
        type_arr.push_back(0.0);
        bound_arr.push_back(x);
        bound_arr.push_back(params.L_d[1]+ params.o_d[1]);
        bound_arr.push_back(z);
        type_arr.push_back(0.0);

        }
    }

    vector<double> o_d = params.o_d, L_d = params.L_d;
    // Shift the center and origin to a get "en quiconce" boundaries
    for (int i = 0; i < 3; i++)
    {
        o_d[i] = params.o_d[i] + params.s * 0.5;
        L_d[i] = params.L_d[i] - params.s;
        // cout << " le centre et longueur de l'axe " << i << "est "<< o_d[i] << " et  " << L_d[i] << endl;
    }

    ni = int(ceil(L_d[0] / params.s));
    dx = params.L_d[0] / ni;
    ++ni;
    nj = int(ceil(L_d[1] / params.s));
    dy = params.L_d[1] / nj;
    ++nj;
    nk= int(ceil(L_d[2] / params.s));
    dz = L_d[2] / nk;
    ++nk;

    // Apply the second layer of FP
    for (int i = 0; i < ni ; ++i) // along x
    {
        double x = o_d[0] + i * dx;
        for (int j = 0; j < nj; ++j) // along y
        {
            double y = params.o_d[1] + j * dy;
            bound_arr.push_back(x);
            bound_arr.push_back(y);
            bound_arr.push_back(params.o_d[2]);
            type_arr.push_back(0.0);
            // bound_arr.push_back(x);
            // bound_arr.push_back(y);
            // bound_arr.push_back(L_d[2] + o_d[2]);
            // type_arr.push_back(0.0);
        }
    }

    for (int j = 0; j < nj; ++j) // along y
    {
        double y = o_d[1] + j * dy;
        for (int k = 1; k < nk-1; ++k) // along z
        {
        
        double z = o_d[2] + k * dz;
        bound_arr.push_back(o_d[0]);
        bound_arr.push_back(y);
        bound_arr.push_back(z);
        type_arr.push_back(0.0);
        bound_arr.push_back(L_d[0] + o_d[0]);
        bound_arr.push_back(y);
        bound_arr.push_back(z);
        type_arr.push_back(0.0);

        }
    }
    for (int i = 1; i < ni-1; ++i) // along x
    {
        double x = o_d[0] + i * dx;
        for (int k = 1; k < nk-1; ++k) // along z
        {
         double z = o_d[2] + k * dz;
        bound_arr.push_back(x);
        bound_arr.push_back(o_d[1]);
        bound_arr.push_back(z);
        type_arr.push_back(0.0);
        bound_arr.push_back(x);
        bound_arr.push_back(L_d[1]+ o_d[1]);
        bound_arr.push_back(z);
        type_arr.push_back(0.0);

        }
    }
}