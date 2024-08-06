#ifndef KERNEL_FUNCTIONS_H
#define KERNEL_FUNCTIONS_H

#include <iostream>
#include <cmath>
#include "structure.h"

double WCoh(double r, GeomData &geomParams, SimulationData simParams);

double WAdh(double r, GeomData &geomParams, SimulationData simParams);

double f_cubic_spline(double r, GeomData &geomParams, SimulationData &simParams);
double derive_cubic_spline(double r, GeomData &geomParams, SimulationData &simParams);

double f_wendland_quintic(double r, GeomData &geomParams, SimulationData &simParams);
double derive_wendland_quintic(double r, GeomData &geomParams, SimulationData &simParams);

#endif // KERNEL_FUNCTIONS_H
