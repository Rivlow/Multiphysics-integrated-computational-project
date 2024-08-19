#include <iostream>
#include <cmath>
#include <stdio.h>
#include "Kernel.h"
#include "structure.h"

using namespace std;

double WCoh(double r, GeomData &geomParams, SimulationData simParams){

    double W = 0.0;
    double coef = simParams.coh_kernel_coef;
    double h = geomParams.h;
    double kappa =1;// geomParams.kappa;
 
    if(r/h <= kappa/2)
        W = coef *(2.0 * (kappa*h-r) * (kappa*h-r) * (kappa*h-r) * r*r*r - kappa*h*kappa*h*kappa*h*kappa*h*kappa*h*kappa*h/64);

    else if(r/h > kappa/2 && r/h <= kappa)
        W = coef*(kappa*h-r)*(kappa*h-r)*(kappa*h-r)*r*r*r;

    else
        W = 0.0;
    
    return W;
}

double WAdh(double r, GeomData &geomParams, SimulationData simParams){

    double W = 0.0;
    double coef = simParams.adh_kernel_coef;
    double h = geomParams.h;
    double kappa = 1;// geomParams.kappa;

    if(r/h > 1/2 && r/h <= kappa)
        W = coef*sqrt(sqrt(-4*r*r/(kappa*h) + 6*r - 2*(kappa*h)));
    else
        W = 0.0;
    
    return W;
}



double CubicSpline(double r, GeomData &geomParams, SimulationData &simParams){

    double alpha = simParams.cubic_kernel_coef;
    double W;
    double h = geomParams.h;
    double kappa = geomParams.kappa;

    if (r/h < kappa/2)
        W = alpha * (2.0/3.0 - r * r / (h * h) + 0.5 * r * r * r / (h * h * h));      
    
    else if (r/h >= kappa/2 && r/h < kappa)
        W = alpha * 1.0/6.0 * ((2.0 - r / h) * (2.0 - r / h) * (2.0 - r / h));
    
    else
        W = 0;
    
    return W;
}

double deriveCubicSpline(double r, GeomData &geomParams, SimulationData &simParams){
    
    double alpha = simParams.cubic_kernel_coef;
    double h = geomParams.h;
    double DW = 0;
    double kappa = geomParams.kappa;

    if (kappa/2 <= r / h && r / h < kappa) 
        DW = alpha / h * (-0.5 * (2.0 - r / h) * (2.0 - r / h));
    
    else if (r / h < kappa/2) 
        DW = alpha / h * (1.5 * r * r / (h * h) - 2.0 * r / h);

    else 
        DW = 0;
        
    return DW;
}

double WendlandQuintic(double r, GeomData &geomParams, SimulationData &simParams){

    double alpha = simParams.quintic_kernel_coef;
    double W = 0.0;
    double h = geomParams.h;
    double kappa = geomParams.kappa;

    if (r/h < kappa) {
        double factor = (1.0 - 0.5 * (r/h));
        W = alpha * (factor*factor*factor*factor) * (2.0 * (r/h) + 1.0);
    }

    return W;
}

double deriveWendlandQuintic(double r, GeomData &geomParams, SimulationData &simParams){

    double alpha = simParams.quintic_kernel_coef;
    double DW = 0.0;
    double h = geomParams.h;
    double kappa = geomParams.kappa;

    if (r/h < kappa) {
        double factor = (1.0 - 0.5* (r/h));
        DW = -5.0*alpha/h * (factor*factor*factor) * (r/h) ;
    }

    return DW;
}
