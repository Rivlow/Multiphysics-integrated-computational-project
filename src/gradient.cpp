#include <stdio.h>
#include <vector>
#include <omp.h>

#include "gradient.h"
#include "find_neighbours.h"
#include "Kernel_functions.h"
#include "tools.h"
#include "structure.h"


using namespace std;

void gradW(const SimulationData& params, 
           vector<vector<double>> &gradW_matrix,
           vector<vector<int>> &neighbours_matrix,
           vector<double> &pos_array){

    double h = params.h;
    int nb_moving_part = params.nb_moving_part;


    // Iterations over each particle
    //#pragma omp parallel for
    for (int pos = 0; pos < int(pos_array.size()/3); pos++)
    {
        vector<int> &neighbours_array = neighbours_matrix[pos];
        // cout << "len(neighbour_list) : " << neighbours_list.size() << endl;

        // Iterations over each associated neighbours of prescribed particles
        for (int idx = 0; idx < int(neighbours_array.size()); idx++)
        {

            int idx_neighbour = neighbours_array[idx];
            double rx, ry, rz, r_ab;
            // cout << "entrance in neighbour loop (before rx, ry, rz) \n";
            rx = (pos_array[3 * pos + 0] - pos_array[3 * idx_neighbour + 0]) * (pos_array[3 * pos + 0] - pos_array[3 * idx_neighbour + 0]);
            ry = (pos_array[3 * pos + 1] - pos_array[3 * idx_neighbour + 1]) * (pos_array[3 * pos + 1] - pos_array[3 * idx_neighbour + 1]);
            rz = (pos_array[3 * pos + 2] - pos_array[3 * idx_neighbour + 2]) * (pos_array[3 * pos + 2] - pos_array[3 * idx_neighbour + 2]);
            r_ab = sqrt(rx + ry + rz);

            // cout << "entrance in neighbour loop (after rx, ry, rz) \n";
            // cout << "Val in neighbour list : (";
            // cout << ")" << endl;
            // cout << "val1 : " << (pos_array[3*pos+0] - pos_array[3*idx_neighbour+0])/r_ab * derive_cubic_spline(r_ab, h)<< endl;
            // cout <<  "val2 : " << (pos_array[3*pos+1] - pos_array[3*idx_neighbour+1])/r_ab * derive_cubic_spline(r_ab, h)<< endl;
            // cout << "val3 : " << (pos_array[3*pos+2] - pos_array[3*idx_neighbour+2])/r_ab * derive_cubic_spline(r_ab, h)<< endl;
            // cout << "pos used : " << pos << " and idx_neighbour used : " << idx_neighbour << endl;
            // cout << "len (pos_array) : " << pos_array.size() << endl;
            // cout << "len (neighbour_list) : " << neighbours_list.size() << endl;
            // cout << "pos_array associated : " << pos_array[3*pos+0] << endl;
            // cout << "idx_neighbour associated : " << pos_array[3*idx_neighbour+0] << endl;

            double val_0 = (pos_array[3 * pos + 0] - pos_array[3 * idx_neighbour + 0]) / r_ab * derive_cubic_spline(r_ab, h);
            double val_1 = (pos_array[3 * pos + 1] - pos_array[3 * idx_neighbour + 1]) / r_ab * derive_cubic_spline(r_ab, h);
            double val_2 = (pos_array[3 * pos + 2] - pos_array[3 * idx_neighbour + 2]) / r_ab * derive_cubic_spline(r_ab, h);

            gradW_matrix[pos].push_back(val_0);
            // cout << "after first push_back"<<endl;
            gradW_matrix[pos].push_back(val_1);
            // cout << "after second push_back"<<endl;
            gradW_matrix[pos].push_back(val_2);
            // cout << "after third push_back"<<endl;
        }
    }

    if (params.PRINT){
            cout << "gradW passed" << endl;
    }
}

void setSpeedOfSound(const SimulationData& params,
                     vector<double> &c_array,
                     vector<double> &rho_array){

    string state_equation = params.state_equation;
    double c_0 = params.c_0;
    double rho_0 = params.rho_0;
    double gamma = params.gamma;

    //#pragma omp parallel for
    for (int pos = 0; pos < int(rho_array.size()); pos++)
    {

        if (state_equation == "Ideal gaz law")
        {
            c_array[pos] = c_0;
        }
        if (state_equation == "Quasi incompresible fluid")
        {
            c_array[pos] = c_0 * pow(rho_array[pos] / rho_0, 0.5 * (gamma - 1));
        }
    }

}

void setPressure(const SimulationData& params,
                 vector<double> &p_array,
                 vector<double> &rho_array){


    string state_equation = params.state_equation;
    double c_0 = params.c_0;
    double rho_0 = params.rho_0;
    double gamma = params.gamma;
    double R = params.R;
    double T = params.T;
    double M = params.M;
    bool PRINT = params.PRINT;
    int nb_moving_part = params.nb_moving_part;

    //#pragma omp parallel for
    for (int pos = 0; pos < nb_moving_part; pos++)
    {

        if (state_equation == "Ideal gaz law")
        {
            p_array[pos] = (rho_array[pos] / rho_0 - 1) * (R * T) / M;
        }

        if (state_equation == "Quasi incompresible fluid")
        {
            double B = c_0 * c_0 * rho_0 / gamma;
            p_array[pos] = B * (pow(rho_array[pos] / rho_0, gamma) - 1);
        }
    }

    if (PRINT){
            cout << "setPressure passed" << endl;
    }
}

void setArtificialViscosity(const SimulationData& params,
                            int t,
                            vector<vector<double>> &artificial_visc_matrix,
                            vector<vector<int>> &neighbours_matrix,
                            vector<double> &c_array,
                            vector<double> &pos_array,
                            vector<double> &rho_array,
                            vector<double> &u_array){

    double beta = params.beta;
    double alpha = params.alpha;
    double h = params.h;
    bool PRINT = params.PRINT;
    int nb_moving_part = params.nb_moving_part;


    if (t == 0)
    {
        //#pragma omp parallel for
        for (int pos = 0; pos < nb_moving_part; pos++)
        {

            vector<int> &neighbours_array = neighbours_matrix[pos];

            for (int idx_neighbour = 0; idx_neighbour < int(neighbours_array.size()); idx_neighbour++)
            {
                artificial_visc_matrix[pos].push_back(0.0);
            }
        }
    }

    else
    {

        vector<double> rel_displ(3), rel_vel(3);

        // Iterations over each particle
        //#pragma omp parallel for
        for (int pos = 0; pos < int(pos_array.size()/3); pos++)
        {

            vector<int> &neighbours_arr = neighbours_matrix[pos];

            // Iteration over each associated neighbours
            for (int idx_neighbour = 0; idx_neighbour < int(neighbours_arr.size()); idx_neighbour++)
            {

                int neighbour_value = neighbours_arr[idx_neighbour];

                rel_displ[0] = (pos_array[3 * pos + 0] - pos_array[3 * neighbour_value + 0]);
                rel_displ[1] = (pos_array[3 * pos + 1] - pos_array[3 * neighbour_value + 1]);
                rel_displ[2] = (pos_array[3 * pos + 2] - pos_array[3 * neighbour_value + 2]);

                rel_vel[0] = (u_array[3 * pos + 0] - u_array[3 * neighbour_value + 0]);
                rel_vel[1] = (u_array[3 * pos + 1] - u_array[3 * neighbour_value + 1]);
                rel_vel[2] = (u_array[3 * pos + 2] - u_array[3 * neighbour_value + 2]);

                double u_ab_x_ab = 0, x_ab_2 = 0;

                // Dot product
                for (int cord = 0; cord < 3; cord++)
                {
                    u_ab_x_ab += rel_vel[cord] * rel_displ[cord];
                    x_ab_2 += rel_displ[cord] * rel_displ[cord];
                }

                double c_a = c_array[pos];
                double c_b = c_array[neighbour_value];
                double rho_a = rho_array[pos];
                double rho_b = rho_array[neighbour_value];

                // cout << "c_a: " << c_a;
                // cout << " c_b: " << c_b << endl;

                double c_ab = 0.5 * (c_a + c_b);
                double rho_ab = 0.5 * (rho_a + rho_b);
                double nu_2 = 0.01 * h * h;
                double mu_ab = (h * u_ab_x_ab) / (x_ab_2 + nu_2);

                // cout << "c_ab: " << c_ab;
                // cout << " rho_ab: " << rho_ab;
                // cout << " nu_2: " << nu_2;
                // cout << " mu_ab: " << mu_ab << "\n " <<endl;

                artificial_visc_matrix[pos].push_back((u_ab_x_ab < 0) ? (-alpha * c_ab * mu_ab + beta * mu_ab * mu_ab) / rho_ab : 0);
            }
        }
    }
    
    if (PRINT){
            cout << "setArtificialViscosity passed" << endl;
    }
}

void continuityEquation(const SimulationData& params,
                        vector<vector<int>> &neighbours_matrix,
                        vector<vector<double>> &gradW_matrix,
                        vector<double> &pos_array,
                        vector<double> &u_array,
                        vector<double> &drhodt_array,
                        vector<double> &rho_array,
                        vector<double> &mass_array){

    bool PRINT = params.PRINT;
    int nb_moving_part = params.nb_moving_part;
             
    // Iterations over each particle
    #pragma omp parallel for
    for (int pos = 0; pos < int(pos_array.size()/3); pos++){

        vector<int> &neighbours_array = neighbours_matrix[pos];
        vector<double> &gradW_array = gradW_matrix[pos];

        // Summation over b = 1 -> nb_neighbours
        for (int idx = 0; idx < int(neighbours_array.size()); idx++){

            // Dot product of u_ab with grad_a(W_ab)
            double dot_product = 0;
            int idx_neighbour = neighbours_array[idx];
            double m_b = mass_array[idx_neighbour];
            for (int cord = 0; cord < 3; cord++){

                double u_a = u_array[3 * pos + cord];
                double u_b = u_array[3 * idx_neighbour + cord];
                dot_product += (u_a - u_b) * gradW_array[3*idx + cord];
            }
            
            drhodt_array[pos] += m_b * dot_product;
        }

    }

    if (PRINT){
            cout << "continuityEquation passed" << endl;
    }
}

void momentumEquation(const SimulationData& params,
                      int t,
                      vector<vector<int>> &neighbours_matrix,
                      vector<vector<double>> &gradW_matrix,
                      vector<vector<double>> &artificial_visc_matrix,
                      vector<double> &mass_array,
                      vector<double> &dudt_array,
                      vector<double> &rho_array,
                      vector<double> &p_array, 
                      vector<double> &c_array,
                      vector<double> &pos_array,
                      vector<double> &u_array){


    double g = params.g;
    bool PRINT = params.PRINT;
    int nb_moving_part = params.nb_moving_part;


    // Compute pressure for all particles
    setPressure(params, p_array, rho_array); 

    // Compute speed of sound for all particles
    setSpeedOfSound(params, c_array, rho_array);

    // Compute artificial viscosity Π_ab for all particles
    setArtificialViscosity(params, t, artificial_visc_matrix, neighbours_matrix,
                            c_array, pos_array, rho_array, u_array); 

    // Iterations over each particle
    #pragma omp parallel for
    for (int pos = 0; pos < nb_moving_part; pos++)
    {
        vector<int> &neighbours_array = neighbours_matrix[pos];
        vector<double> &gradW_array = gradW_matrix[pos];
        vector<double> &artificial_visc_array = artificial_visc_matrix[pos];
        vector<double> F_vol = {0.0, 0.0, g};
        double rho_a = rho_array[pos];
        double p_a = p_array[pos];

        for (int cord = 0; cord < 3; cord++)
        {
            // Summation over b = 1 -> nb_neighbours
            for (int idx_neighbour = 0; idx_neighbour < int(neighbours_array.size()); idx_neighbour++)
            {
                double pi_ab = artificial_visc_array[idx_neighbour];
                double rho_b = rho_array[neighbours_array[idx_neighbour]];
                double m_b = mass_array[neighbours_array[idx_neighbour]];
                double p_b = p_array[neighbours_array[idx_neighbour]];

                dudt_array[3 * pos + cord] += m_b * (p_b / (rho_b * rho_b) + p_a / (rho_a * rho_a) + pi_ab) 
                                            * gradW_array[3*idx_neighbour + cord];

            }

            dudt_array[3 * pos + cord] *= -1;
            dudt_array[3 * pos + cord] += F_vol[cord];
        }
    }

    if (PRINT){
        cout << "momentumEquation passed" << endl;
    }
}

void update(const SimulationData& params,
            vector<double> &pos_array,
            vector<double> &u_array,
            vector<double> &rho_array,
            vector<double> &drhodt_array,
            vector<double> &dudt_array){

    double dt = params.dt;
    bool PRINT = params.PRINT;
    
    for (int pos = 0; pos < int(pos_array.size()/3); pos++){

        rho_array[pos] += dt * drhodt_array[pos];

        for (int cord = 0; cord < 3; cord++){

            pos_array[3 * pos + cord] += dt * u_array[3 * pos + cord];
            u_array[3 * pos + cord] += dt * dudt_array[3 * pos + cord];
        }
    }

    if (PRINT){
        cout << "update passed" << endl;
    }
}
