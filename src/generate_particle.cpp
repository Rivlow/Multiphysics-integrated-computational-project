#include "generate_particle.h"
#include <iostream>
#include <cmath> // ceil
using namespace std;

/**
 * @brief build a cube of particles aligned with x,y,z axes.
 *
 * @param o corner of the cube with the lowest (x,y,z) values
 * @param L edge lengths along x,y and z
 * @param s particle spacing
 * @param pos_arr pos_arritions of particles
 */

int evaluateNumberParticles(vector<double> &L, 
                            double &s){

    int ni = int(ceil(L[0] / s));
    // double dx = L[0] / ni;
    ++ni;
    int nj = int(ceil(L[1] / s));
    // double dy = L[1] / nj;
    ++nj;
    int nk = int(ceil(L[2] / s));
    // double dz = L[2] / nk;
    ++nk;

    std::cout << "\t=> " << ni << "*" << nj << "*" << nk << " = " << ni * nj * nk << " particles to be generated\n";

    return ni * nj * nk;
}

void meshcube(vector<double> &o, 
              vector<double> &L, 
              vector<double> &pos_arr,
              vector<double> &type_arr,
              double s){

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

                type_arr.push_back(1.0);
            }
        }
    }
}

void meshBoundary(vector<double> &o_d, 
                  vector<double> &L_d, 
                  vector<double> &bound_arr, 
                  vector<double> &type_arr,
                  double s){

    int ni = int(ceil(L_d[0] / s));
    double dx = L_d[0] / ni;
    ++ni;
    int nj = int(ceil(L_d[1] / s));
    double dy = L_d[1] / nj;
    ++nj;
    int nk = int(ceil(L_d[2] / s));
    double dz = L_d[2] / nk;
    ++nk;
    cout << "ni, nj, nk = " << ni << nj << nk << endl;

    bound_arr.reserve(bound_arr.size() + 6 * (ni * nj + (nk - 2) * nj + (ni - 2) * (nk - 2)));

    // Apply first layer of FP
    for (int i = 0; i < ni ; ++i) // along x
    {
        double x = o_d[0] + i * dx;
        for (int j = 0 ; j < nj ; ++j) // along y
        {
            double y = o_d[1] + j * dy;
            bound_arr.push_back(x);
            bound_arr.push_back(y);
            bound_arr.push_back(o_d[2]);
            type_arr.push_back(0.0);
            bound_arr.push_back(x);
            bound_arr.push_back(y);
            bound_arr.push_back(L_d[2] + o_d[2]);
            type_arr.push_back(0.0);
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

    // Shift the center and origin to a get "en quiconce" boundaries
    for (size_t i = 0; i < 3; i++)
    {
        o_d[i] = o_d[i] + s * 0.5;
        L_d[i] = L_d[i] - s;
        // cout << " le centre et longueur de l'axe " << i << "est "<< o_d[i] << " et  " << L_d[i] << endl;
    }

    ni = int(ceil(L_d[0] / s));
    dx = L_d[0] / ni;
    ++ni;
    nj = int(ceil(L_d[1] / s));
    dy = L_d[1] / nj;
    ++nj;
    nk= int(ceil(L_d[2] / s));
    dz = L_d[2] / nk;
    ++nk;

    // Apply the second layer of FP
    for (int i = 0; i < ni ; ++i) // along x
    {
        double x = o_d[0] + i * dx;
        for (int j = 0; j < nj; ++j) // along y
        {
            double y = o_d[1] + j * dy;
            bound_arr.push_back(x);
            bound_arr.push_back(y);
            bound_arr.push_back(o_d[2]);
            type_arr.push_back(0.0);
            bound_arr.push_back(x);
            bound_arr.push_back(y);
            bound_arr.push_back(L_d[2] + o_d[2]);
            type_arr.push_back(0.0);
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
   /* for (size_t i = 0; i < 3; i++)
    {
        o_d[i] = o_d[i] + s * 0.5;
        L_d[i] = L_d[i] - s;
        // cout << " le centre et longueur de l'axe " << i << "est "<< o_d[i] << " et  " << L_d[i] << endl;
    }
    
    ni = int(ceil(L_d[0] / s));
    dx = L_d[0] / ni;
    ++ni;
    nj = int(ceil(L_d[1] / s));
    dy = L_d[1] / nj;
    ++nj;
    nk= int(ceil(L_d[2] / s));
    dz = L_d[2] / nk;
    ++nk;

    // Apply the second layer of FP
    for (int i = 0; i < ni ; ++i) // along x
    {
        double x = o_d[0] + i * dx;
        for (int j = 0; j < nj; ++j) // along y
        {
            double y = o_d[1] + j * dy;
            bound_arr.push_back(x);
            bound_arr.push_back(y);
            bound_arr.push_back(o_d[2]);
            type_arr.push_back(0.0);
            bound_arr.push_back(x);
            bound_arr.push_back(y);
            bound_arr.push_back(L_d[2] + o_d[2]);
            type_arr.push_back(0.0);
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
    }*/
}