#ifndef KERNEL_FUNCTIONS_H
#define KERNEL_FUNCTIONS_H

#include <iostream>
#include <cmath>
#include "structure.h"

double WCoh(double r, GeomData &geomParams, SimulationData simParams);

double WAdh(double r, GeomData &geomParams, SimulationData simParams);

double CubicSpline(double r, GeomData &geomParams, SimulationData &simParams);
double deriveCubicSpline(double r, GeomData &geomParams, SimulationData &simParams);

double WendlandQuintic(double r, GeomData &geomParams, SimulationData &simParams);
double deriveWendlandQuintic(double r, GeomData &geomParams, SimulationData &simParams);

#endif // KERNEL_FUNCTIONS_H
