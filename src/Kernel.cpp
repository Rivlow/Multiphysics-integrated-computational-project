#include <iostream>
#include <cmath>
#include <stdio.h>
#include "Kernel.h"
#include "structure.h"

using namespace std;

double W_coh(double r, GeomData &geomParams, SimulationData simParams){

    double W = 0.0;
    double coef = simParams.coh_kernel_coef;
    double h = geomParams.h;

    if(2.0*r>h && r<=h)
        W = coef*(h-r)*(h-r)*(h-r)*r*r*r;
    
    else if(r>0 && 2.0*r<=h)
        W = coef *( 2.0*(h-r)*(h-r)*(h-r)*r*r*r - h*h*h*h*h*h/64);
    
    else
        W = 0.0;
    
    return W;
}

double W_adh(double r, GeomData &geomParams, SimulationData simParams){

    double W = 0.0;
    double coef = simParams.adh_kernel_coef;
    double kh = geomParams.kappa*geomParams.h;;

    if(2.0*r>kh && r<=kh)
        W = coef*sqrt(sqrt(-4*r*r/kh+6*r-2*kh));
    else
        W = 0.0;
    
    return W;
}



double f_cubic_spline(double r, GeomData &geomParams, SimulationData &simParams){

    double alpha = simParams.cubic_kernel_coef;
    double W;
    double h = geomParams.h;

    if (r / h < 1.0)
        W = alpha * (2.0/3.0 - r * r / (h * h) + 0.5 * r * r * r / (h * h * h));      
    
    else if (1.0 <= r / h && r / h < 2.0)
        W = alpha * 1.0/6.0 * ((2.0 - r / h) * (2.0 - r / h) * (2.0 - r / h));
    
    else
        W = 0;
    
    return W;
}

double derive_cubic_spline(double r, GeomData &geomParams, SimulationData &simParams){
    
    double alpha = simParams.cubic_kernel_coef;
    double h = geomParams.h;
    double DW = 0;

    if (1.0 <= r / h && r / h < 2.0) 
        DW = alpha / h * (-0.5 * (2.0 - r / h) * (2.0 - r / h));
    
    else if (r / h < 1.0) 
        DW = alpha / h * (1.5 * r * r / (h * h) - 2.0 * r / h);

    else 
        DW = 0;
        
    return DW;
}
