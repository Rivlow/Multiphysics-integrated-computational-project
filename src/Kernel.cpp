#include <iostream>
#include <cmath>
#include <stdio.h>
#include "Kernel.h"
#include "structure.h"

using namespace std;

double W_coh(double r, double h, SimulationData simParams){
    double W = 0.0;
    double cst = (simParams.dimension == 3) ? 32/(M_PI*h*h*h*h*h*h*h*h*h) : 40/(M_PI*h*h*h*h*h*h*h*h);
    if(2.0*r>h && r<=h){
        W = cst*(h-r)*(h-r)*(h-r)*r*r*r;
    }
    else if(r>0 && 2.0*r<=h){
        W = cst *( 2.0*(h-r)*(h-r)*(h-r)*r*r*r - h*h*h*h*h*h/64);
    }
    else{
        W = 0.0;
    }
    return W;
}

double W_adh(double r, double h, SimulationData simParams){

    double W = 0.0;
    double cst = (simParams.dimension == 2)? 16/(4* M_PI*pow(h, 2.25)): 0.0007/pow(h,3.25);

    if(2.0*r>h && r<=h)
        W = cst*sqrt(sqrt(-4*r*r/h+6*r-2*h));
    else
        W = 0.0;
    
    return W;
}



double f_gaussian(double r, double h)
{
    double alpha = 1 / (M_PI * sqrt(M_PI) * h * h * h);
    double W = alpha * exp(-r * r / (h * h));
    return W;
}

double deriv_gaussian(double r, double h){

    double alpha = 1 / (M_PI * sqrt(M_PI) * h * h * h);
    double DW = alpha / h * (-2 * r / h * exp(-r * r / (h * h)));
    return DW;
}

double f_bell(double r, double h){

    double alpha = 105 / (16 * M_PI * h * h * h);
    double W;

    if (r / h <= 1)
        W = alpha * (1 + 3 * r / h) * (1 - r / h) * (1 - r / h) * (1 - r / h);
    else
        W = 0;
    
    return W;
}

double derive_bell(double r, double h){

    double alpha = 105 / (16 * M_PI * h * h * h);
    double DW;

    if (r / h <= 1)
        DW = alpha / h * 3 * ((1 - r / h) * (1 - r / h) * (1 - r / h) - (1 + 3 * r / h) * (1 - r / h) * (1 - r / h) * (1 - r / h));
    else
        DW = 0;
    return DW;
}

double f_cubic_spline(double r, double h, SimulationData &simParams){
    int dim = simParams.dimension;
    double alpha = 0;
    if(dim == 3){
        alpha = 3.0 / (2.0 * M_PI * h * h * h);
    }
    if(dim == 2){
        alpha = 15.0 /( 7.0 * M_PI * h * h );
    }
    double W;

    if (r / h < 1.0){
        W = alpha * (2.0/3.0 - r * r / (h * h) + 0.5 * r * r * r / (h * h * h));
        
    }
        
    
    else if (1.0 <= r / h && r / h < 2.0)
        W = alpha * 1.0/6.0 * ((2.0 - r / h) * (2.0 - r / h) * (2.0 - r / h));
    
    else
        W = 0;
    
    return W;
}

double derive_cubic_spline(double r, double h, SimulationData &simParams)
{

    int dim = simParams.dimension;
    double alpha = 0;
    if(dim == 3){
        alpha = 3.0 / (2.0 * M_PI * h * h * h);
    }
    if(dim == 2){
        alpha = 15.0 /( 7.0 * M_PI * h * h );
    }
    
    double DW = 0;

    if (1.0 <= r / h && r / h < 2.0) 
        DW = alpha / h * (-0.5 * (2.0 - r / h) * (2.0 - r / h));
    
    else if (r / h < 1.0) 
        DW = alpha / h * (1.5 * r * r / (h * h) - 2.0 * r / h);

    else 
        DW = 0;
        
    return DW;
}

double f_quadratic(double r, double h){

    double alpha = 5 / (4 * M_PI * h * h * h);
    double W = 0;

    if (r / h <= 2)
        W = alpha * (3 / 16 * r * r / (h * h) - 3 / 4 * r / h + 3 / 4);
    
    return W;
}

double derive_quadratic(double r, double h){

    double alpha = 5 / (4 * M_PI * h * h * h);
    double DW;

    if (r / h <= 2)
        DW = alpha / h * (3 / 8 * r / h - 3 / 4);
    else
        DW = 0;
    
    return DW;
}

double f_quintic(double r, double h){

    double alpha = 21 / (16 * M_PI * h * h * h);
    double W;

    if (r / h <= 2)
        W = alpha * ((1 - 1 / 2 * r / h) * (1 - 1 / 2 * r / h) * (1 - 1 / 2 * r / h) * (1 - 1 / 2 * r / h) * (2 * r / h + 1));
    else
        W = 0;
    
    return W;
}

double derive_quintic(double r, double h){

    double alpha = 21 / (16 * M_PI * h * h * h);
    double DW;

    if (r / h <= 2)
        DW = alpha / h * (2 * (1 - 0.5 * r / h) * (1 - 0.5 * r / h) * (1 - 0.5 * r / h) * (1 - 0.5 * r / h) - (2 * r / h - 1) * (1 - 0.5 * r / h) * (1 - 0.5 * r / h) * (1 - 0.5 * r / h));
    else
        DW = 0;

    return DW;
}

double f_quinitc_spline(double r, double h){

    double alpha = 3 / (359 * M_PI * h * h * h);
    double W = 0;
    
    if (r / h < 1)
        W = alpha * ((3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) - 6 * (2 - r / h) * (2 - r / h) * (2 - r / h) * (2 - r / h) * (2 - r / h) + 15 * (1 - r / h) * (1 - r / h) * (1 - r / h) * (1 - r / h) * (1 - r / h));
    
    if (1 <= r / h && r / h < 2)
        W = alpha * ((3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) - 6 * (2 - r / h) * (2 - r / h) * (2 - r / h) * (2 - r / h) * (2 - r / h));
    
    if (2 <= r / h && r / h < 3)
        W = alpha * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h);
    
    return W;
}

double derive_quintic_spline(double r, double h){

    double alpha = 3 / (359 * M_PI * h * h * h);
    double DW;

    if (r / h < 1)
        DW = alpha / h * (-5 * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) + 30 * (2 - r / h) * (2 - r / h) * (2 - r / h) * (2 - r / h) - 75 * (1 - r / h) * (1 - r / h) * (1 - r / h) * (1 - r / h));
    
    else if (1 <= r / h && r / h < 2)
        DW = alpha / h * (-5 * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) + 30 * (2 - r / h) * (2 - r / h) * (2 - r / h) * (2 - r / h));
    
    else if (2 <= r / h && r / h < 3)
        DW = -5 * alpha / h * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h);

    else
        DW = 0;
    
    return DW;
}
