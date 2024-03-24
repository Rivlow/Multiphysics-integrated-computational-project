#include "gradient.h"
#include <stdio.h>
#include <vector>
#include "sorted_list.h"
#include "Kernel_functions.h"
using namespace std;



void gradW(vector<vector<double>> &gradW_matrix, int &nb_moving_part, vector<double> &pos_arr, vector<vector<int>> &neighbours_matrix, 
                    double &h, int &Nx, int &Ny, int &Nz){
    
    // Iterations over each particle
    for (size_t pos = 0; pos < nb_moving_part; pos++){
        vector<int> &neighbours_list = neighbours_matrix[pos];
        //cout << "len(neighbour_list) : " << neighbours_list.size() << endl;

        // Iterations over each associated neighbours of prescribed particles
        for (size_t idx = 0; idx < neighbours_list.size(); idx++){     
            
            size_t idx_neighbour = neighbours_list[idx];
            double rx, ry, rz, r_ab;
            //cout << "entrance in neighbour loop (before rx, ry, rz) \n";
            rx = (pos_arr[3*pos+0] - pos_arr[3*idx_neighbour+0])*(pos_arr[3*pos+0] - pos_arr[3*idx_neighbour+0]);
            ry = (pos_arr[3*pos+1] - pos_arr[3*idx_neighbour+1])*(pos_arr[3*pos+1] - pos_arr[3*idx_neighbour+1]);
            rz = (pos_arr[3*pos+2] - pos_arr[3*idx_neighbour+2])*(pos_arr[3*pos+2] - pos_arr[3*idx_neighbour+2]);
            r_ab = sqrt(rx + ry + rz);

            //cout << "entrance in neighbour loop (after rx, ry, rz) \n";
            //cout << "Val in neighbour list : (";
            for (size_t idxx = 0; idxx < neighbours_list.size(); idxx++){
                //cout << neighbours_list[idxx] << " , ";
            }

            //cout << ")" << endl;
            //cout << "val1 : " << (pos_arr[3*pos+0] - pos_arr[3*idx_neighbour+0])/r_ab * derive_cubic_spline(r_ab, h)<< endl;
            //cout <<  "val2 : " << (pos_arr[3*pos+1] - pos_arr[3*idx_neighbour+1])/r_ab * derive_cubic_spline(r_ab, h)<< endl;
            //cout << "val3 : " << (pos_arr[3*pos+2] - pos_arr[3*idx_neighbour+2])/r_ab * derive_cubic_spline(r_ab, h)<< endl;
            //cout << "pos used : " << pos << " and idx_neighbour used : " << idx_neighbour << endl;
            //cout << "len (pos_arr) : " << pos_arr.size() << endl;
            //cout << "len (neighbour_list) : " << neighbours_list.size() << endl;
            //cout << "pos_arr associated : " << pos_arr[3*pos+0] << endl;
            //cout << "idx_neighbour associated : " << pos_arr[3*idx_neighbour+0] << endl;
            
            double val_0 = abs(pos_arr[3*pos+0] - pos_arr[3*idx_neighbour+0])/r_ab * derive_cubic_spline(r_ab, h);
            double val_1 = abs(pos_arr[3*pos+1] - pos_arr[3*idx_neighbour+1])/r_ab * derive_cubic_spline(r_ab, h);
            double val_2 = abs(pos_arr[3*pos+2] - pos_arr[3*idx_neighbour+2])/r_ab * derive_cubic_spline(r_ab, h);
            //cout << "derive cubic spline : " << derive_cubic_spline(r_ab, h) << " r_ab/h : " << r_ab/h << endl;
            gradW_matrix[pos].push_back(val_0);
            //cout << "after first push_back"<<endl;
            gradW_matrix[pos].push_back(val_1);
            //cout << "after second push_back"<<endl;
            gradW_matrix[pos].push_back(val_2);
            //cout << "after third push_back"<<endl;

        }
    }
}


double setSpeedOfSound(double &rho, double &rho_0, double &c_0, double &gamma, string &state_equation_chosen){

    double c = 0;
    
    if (state_equation_chosen == "Ideal gaz law"){
        c = c_0;
    }
    if (state_equation_chosen == "Quasi incompresible fluid"){
        c = c_0*pow(rho/rho_0, 0.5*(gamma-1));
    }
    return c;
}

void setPressure(int &nb_moving_part, vector<double> &p_arr, vector<double> &rho_arr, double &rho_0,  double &c_0, double &R,  double &T,
                      double &M, double &gamma, string &state_equation_chosen){

    for (size_t pos = 0; pos < p_arr.size(); pos++){

        if (state_equation_chosen == "Ideal gaz law"){
            p_arr[pos] = (rho_arr[pos]/rho_0 - 1)*(R*T)/M;
        }

        if (state_equation_chosen == "Quasi incompresible fluid"){
            double B = c_0*c_0*rho_0/gamma;
            p_arr[pos] = B*(pow(rho_arr[pos]/rho_0, gamma) - 1);
        }
    }
}

void setArtificialViscosity(int &t, vector<vector<double>> &artificial_visc, int &nb_moving_part, vector<double> &pos_arr, vector<vector<int>> &neighbours_matrix, vector<double> &rho_arr, 
                            vector<double> &u_arr, double &alpha, double &beta, double &rho_0, double &c_0, double &gamma, double &R, double &T, double &M, double &h, string &state_equation_chosen){

    if (t == 0){
        for (size_t pos = 0; pos < nb_moving_part; pos++){

            vector<int> &neighbours_list = neighbours_matrix[pos];
        
            for (size_t idx_neighbour = 0; idx_neighbour < neighbours_list.size(); idx_neighbour++){
                artificial_visc[pos].push_back(0.0);
            }
        }
    }

    
    else{

        vector<double> rel_displ(3), rel_vel(3);

        // Iterations over each particle
        for (size_t pos = 0; pos < nb_moving_part; pos++){

            vector<int> &neighbours_arr = neighbours_matrix[pos];
           // cout << "For particle : " << pos << endl;
           // cout << " pos (x): " << pos_arr[3*pos+0];
           // cout << " pos (y): " << pos_arr[3*pos+1];
           // cout << " pos (z): " << pos_arr[3*pos+2] << "\n" <<endl;
            
            // Iteration over each associated neighbours
            for (size_t idx_neighbour = 0; idx_neighbour < neighbours_arr.size(); idx_neighbour++){   

                int neighbour_value = neighbours_arr[idx_neighbour];

                rel_displ[0] = (pos_arr[3*pos+0] - pos_arr[3*neighbour_value+0]);
                rel_displ[1] = (pos_arr[3*pos+1] - pos_arr[3*neighbour_value+1]);
                rel_displ[2] = (pos_arr[3*pos+2] - pos_arr[3*neighbour_value+2]);

                rel_vel[0] = (u_arr[3*pos+0] - u_arr[3*neighbour_value+0]);
                rel_vel[1] = (u_arr[3*pos+1] - u_arr[3*neighbour_value+1]);
                rel_vel[2] = (u_arr[3*pos+2] - u_arr[3*neighbour_value+2]);

                
                

                //cout << "vel nei (x): " << u_arr[3*neighbour_value+0];
                //cout << " vel nei (y): " << u_arr[3*neighbour_value+1];
                //cout << " vel nei (z): " << u_arr[3*neighbour_value+2] << endl;
                

                double res = 0, xa_xb = 0;

                // Dot product
                for (size_t cord = 0; cord < 3; cord++){
                    res += rel_vel[cord]*rel_displ[cord];
                    xa_xb += rel_displ[cord]*rel_displ[cord];
                }

                //cout << "res (y): " << res;
                //cout << "xa_xb (z): " << xa_xb << endl;

                double c_a, c_b;
                c_a = setSpeedOfSound(rho_arr[pos], rho_0, c_0, gamma, state_equation_chosen);
                c_b = setSpeedOfSound(rho_arr[neighbour_value], rho_0, c_0, gamma, state_equation_chosen);

                //cout << "c_a: " << c_a;
                //cout << " c_b: " << c_b << endl;

                double c_ab = 0.5*(c_a + c_b);
                double rho_ab = 0.5*(rho_arr[pos] + rho_arr[neighbour_value]);
                double nu_2 = 0.01*h*h;
                double mu_ab = (h*res)/(xa_xb + nu_2);

                //cout << "c_ab: " << c_ab;
                //cout << " rho_ab: " << rho_ab;
                //cout << " nu_2: " << nu_2;
                //cout << " mu_ab: " << mu_ab << "\n " <<endl;

                
                artificial_visc[pos].push_back((res < 0) ? (-alpha*c_ab*mu_ab + beta*mu_ab*mu_ab)/rho_ab : 0);
            }
        }
    }
}

void continuityEquation(int &nb_moving_part, vector<double> &pos_arr, vector<double> &u_arr, vector<vector<int>> &neighbours_matrix, 
                         vector<vector<double>> &gradW_matrix, vector<double> &drhodt_arr, vector<double> &rho_arr, vector<double> &mass_arr,  double &h){

    // Iterations over each particle                    
    for (size_t pos = 0; pos < nb_moving_part; pos++){

        drhodt_arr[pos] = 0;

        vector<int> neighbours_list = neighbours_matrix[pos];
        vector<double> gradW_list = gradW_matrix[pos];
        
        //cout << "for part : " << pos << "number of nei : " << gradW_list.size()  << endl;
        //cout << "for part : " << pos << "number of nei : " << neighbours_list.size()  << endl;
        double drhodt = 0;
        
        // Summation over b = 1 -> nb_neighbours
        for (size_t idx_neighbour = 0; idx_neighbour < neighbours_list.size(); idx_neighbour++){   

            // Dot product of u_ab with grad_a(W_ab)
            double dot_product = 0;
            for (size_t x = 0; x < 3; x++){

                dot_product += (u_arr[3*pos+x] - u_arr[3*neighbours_list[idx_neighbour]+x])*(gradW_list[3*idx_neighbour+x]);
                
            }
            if(u_arr[3*pos+2] - u_arr[3*neighbours_list[idx_neighbour]+2]){
               double rx = (pos_arr[3*pos+0] - pos_arr[3*neighbours_list[idx_neighbour]+0])*(pos_arr[3*pos+0] - pos_arr[3*neighbours_list[idx_neighbour]+0]);
            double ry = (pos_arr[3*pos+1] - pos_arr[3*neighbours_list[idx_neighbour]+1])*(pos_arr[3*pos+1] - pos_arr[3*neighbours_list[idx_neighbour]+1]);
            double rz = (pos_arr[3*pos+2] - pos_arr[3*neighbours_list[idx_neighbour]+2])*(pos_arr[3*pos+2] - pos_arr[3*neighbours_list[idx_neighbour]+2]);
            double r_ab = sqrt(rx + ry + rz);
            //cout << "relative velocity : " << u_arr[3*pos+2] - u_arr[3*neighbours_list[idx_neighbour]+2]<< "relative position : "<< r_ab/h << " and gradW : " << gradW_list[idx_neighbour+2] << endl;

            } 
            drhodt += mass_arr[neighbours_list[idx_neighbour]]*dot_product;
            //if(u_arr[3*pos+2] - u_arr[3*neighbours_list[idx_neighbour]+2]) cout << "drhodt " << drhodt << endl;
            //cout << "mass used : " << mass_arr[neighbours_list[idx_neighbour]] << " and dot product : " << dot_product << endl;
        }

        drhodt_arr[pos] = drhodt;
        
        //cout << "drhost : " << drhodt << endl;
       
    }
}

void momentumEquation(int &nb_moving_part, vector<vector<int>> &neighbours_matrix, vector<double> &mass_arr, vector<vector<double>> &gradW_matrix, 
                      vector<double> &dudt_arr, vector<vector<double>> &artificial_visc_matrix, vector<double> &rho_arr, double &rho_0, double &c_0,
                      vector<double> &p_arr, double &R, double &T, double &M, double &gamma, double &g, string &state_equation_chosen){

    // Iterations over each particle
 
    for (size_t pos = 0; pos < nb_moving_part; pos++){
        
        for (size_t cord = 0; cord < 3; cord++){
        dudt_arr[3*pos+cord] = 0.0;
        }   


        vector<int> &neighbours_arr = neighbours_matrix[pos];
        vector<double> &gradW_arr = gradW_matrix[pos];
        vector<double> &artificial_visc_arr = artificial_visc_matrix[pos];
        vector<double> dudt(3, 0.0), F_vol = {0.0, 0.0, g};



        for (size_t cord = 0; cord < 3; cord++){

            // Summation over b = 1 -> nb_neighbours
            for (size_t idx_neighbour = 0; idx_neighbour < neighbours_arr.size(); idx_neighbour++){   

                double rho_a = rho_arr[pos], rho_b = rho_arr[neighbours_arr[idx_neighbour]];
                double pi_ab = 0; //artificial_visc_arr[idx_neighbour];
                double m_b = mass_arr[neighbours_arr[idx_neighbour]];
                double p_a = p_arr[pos], p_b = p_arr[neighbours_arr[idx_neighbour]];

                dudt[cord] += m_b*(p_b/(rho_b*rho_b) + p_a/(rho_a*rho_a) + pi_ab)*gradW_arr[3*idx_neighbour+cord];/*
                cout << "p_a" << p_a << endl;
                cout << "p_b" << p_b << endl;
                cout << "rho_b" << rho_b << endl;
                cout << "rho_a" << rho_a << endl;
                cout << "m_b" << m_b << endl;*/

            }

            dudt[cord] +=  F_vol[cord];

        }

        for (size_t cord = 0; cord < 3; cord++) {
            dudt[cord] *= -1; 
            //cout << "dudt : " << dudt[cord] << " at cord : " << cord << endl;
            dudt_arr[3*pos+cord] = dudt[cord];
            //cout << "dudt_arr[3*pos+cord] " << dudt_arr[3*pos+cord] << endl;
        }

    }
}

