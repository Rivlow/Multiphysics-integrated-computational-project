#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
#include <cstdlib>
#include <random>
#include <list>
#include <unordered_map>
#include <iostream>
#include <chrono>
#include <algorithm>

#include "sorted_list.h"
#include "integration.h"
#include <omp.h>

void Euler(const double dt, 
           const size_t nb_tot_part,
           vector<double> &rho_array,
           vector<double> &pos_array,
           vector<double> &u_array,
           vector<double> drhodt_array,
           vector<double> dudt_array){
    #pragma omp parallel for   
    for (size_t pos = 0; pos < nb_tot_part; pos++){

            rho_array[pos] += dt * drhodt_array[pos];

            for (size_t cord = 0; cord < 3; cord++)
            {

                pos_array[3 * pos + cord] += dt * u_array[3 * pos + cord];
                u_array[3 * pos + cord] += dt * dudt_array[3 * pos + cord];
            }
        }
}

void runchKutta(const double dt, 
           const size_t nb_tot_part,
           const double theta,
           vector<double> &rho_array,
           vector<double> &pos_array,
           vector<double> &u_array,
           vector<double> drhodt_array,
           vector<double> drhodt_array_half,
           vector<double> dudt_array,
           vector<double> dudt_array_half){

        #pragma omp parallel for
        for (size_t pos = 0; pos < nb_tot_part; pos++){

            rho_array[pos] += dt * ((1-theta)*drhodt_array[pos] + theta*drhodt_array_half[pos]);

            for (size_t cord = 0; cord < 3; cord++)
            {
                double u_temp = u_array[3 * pos + cord] + dt*dudt_array_half[3 * pos + cord];
                pos_array[3 * pos + cord] += dt * ((1-theta)*u_array[3 * pos + cord] + theta*u_temp);
                u_array[3 * pos + cord] += dt *  ((1-theta)*dudt_array[3 * pos + cord] + theta*dudt_array_half[3 * pos + cord]);

            }
        }    
}