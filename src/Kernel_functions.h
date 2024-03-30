#ifndef KERNEL_FUNCTIONS_H
#define KERNEL_FUNCTIONS_H

#include <iostream>
#include <cmath>

double f_gaussian(double r, double h);
double deriv_gaussian(double r, double h);

double f_bell(double r, double h);
double derive_bell(double r, double h);

double f_cubic_spline(double r, double h);
double derive_cubic_spline(double r, double h);

double f_quadratic(double r, double h);
double derive_quadratic(double r, double h);

double f_quintic(double r, double h);
double derive_quintic(double r, double h);

double f_quinitc_spline(double r, double h);
double derive_quintic_spline(double r, double h);

#endif // KERNEL_FUNCTIONS_H
