#include "generate_particle.h"
#include <iostream>
#include <cmath> // ceil

/**
 * @brief build a cube of particles aligned with x,y,z axes.
 *
 * @param o corner of the cube with the lowest (x,y,z) values
 * @param L edge lengths along x,y and z
 * @param s particle spacing
 * @param pos positions of particles
 */

void meshcube(double o[3], double L[3], double s, std::vector<double> &particle_x, 
                std::vector<double> &particle_y, std::vector<double> &particle_z)
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
    std::cout << "\t=> " << ni << "*" << nj << "*" << nk << " = " << ni * nj * nk << " particles to be generated\n";

    /* memory allocation
    particle_x.resize(particle_x.size() + ni * nj * nk);
    particle_y.resize(particle_y.size() + ni * nj * nk);
    particle_z.resize(particle_z.size() + ni * nj * nk);*/

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
                particle_x.push_back(x);
                particle_y.push_back(y);
                particle_z.push_back(z);
            }
        }
    }
}
