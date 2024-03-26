#include <iostream>
#include <cmath>
#include <stdio.h>
#include "Kernel_functions.h"

using namespace std;

double f_gaussian(double r, double h)
{
    double alpha = 1 / (M_PI * sqrt(M_PI) * h * h * h);
    double W = alpha * exp(-r * r / (h * h));
    return W;
}

double deriv_gaussian(double r, double h)
{
    double alpha = 1 / (M_PI * sqrt(M_PI) * h * h * h);
    double DW = alpha / h * (-2 * r / h * exp(-r * r / (h * h)));
    return DW;
}

double f_bell(double r, double h)
{
    double alpha = 105 / (16 * M_PI * h * h * h);
    double W = 0;
    if (r / h <= 1)
    {
        W = alpha * (1 + 3 * r / h) * (1 - r / h) * (1 - r / h) * (1 - r / h);
    }
    return W;
}

double derive_bell(double r, double h)
{
    double alpha = 105 / (16 * M_PI * h * h * h);
    double DW = 0;
    if (r / h <= 1)
    {
        DW = alpha / h * 3 * ((1 - r / h) * (1 - r / h) * (1 - r / h) - (1 + 3 * r / h) * (1 - r / h) * (1 - r / h) * (1 - r / h));
    }
    return DW;
}

double f_cubic_spline(double r, double h)
{
    double alpha = 3 / (2 * M_PI * h * h * h);
    double W = 0;
    if (r / h < 1)
    {
        W = alpha * (3 / 2 - r * r / (h * h) + 1 / 2 * r * r * r / (h * h * h));
    }
    if (1 <= r / h && r / h < 2)
    {

        W = alpha * 1 / 6 * ((1 - r / h) * (1 - r / h) * (1 - r / h));
    }
    return W;
}

double derive_cubic_spline(double r, double h)
{

    double alpha = 3 / (2 * M_PI * h * h * h);
    double DW = 0;

    

    if (1.0 <= r / h && r / h < 2.0)
    {
       DW = alpha / h * (-0.5 * (1 - r / h) * (1 - r / h));
    }
    if (r / h < 1.0)
    {
         DW = alpha / h * (1.5 * r * r / (h * h) - 2 * r / h);
    }

    return DW;
}

double f_quadratic(double r, double h)
{
    double alpha = 5 / (4 * M_PI * h * h * h);
    double W = 0;
    if (r / h <= 2)
    {
        W = alpha * (3 / 16 * r * r / (h * h) - 3 / 4 * r / h + 3 / 4);
    }
    return W;
}

double derive_quadratic(double r, double h)
{
    double alpha = 5 / (4 * M_PI * h * h * h);
    double DW = 0;
    if (r / h <= 2)
    {
        DW = alpha / h * (3 / 8 * r / h - 3 / 4);
    }
    return DW;
}

double f_quintic(double r, double h)
{
    double alpha = 21 / (16 * M_PI * h * h * h);
    double W = 0;
    if (r / h <= 2)
    {
        W = alpha * ((1 - 1 / 2 * r / h) * (1 - 1 / 2 * r / h) * (1 - 1 / 2 * r / h) * (1 - 1 / 2 * r / h) * (2 * r / h + 1));
    }

    return W;
}

double derive_quintic(double r, double h)
{
    double alpha = 21 / (16 * M_PI * h * h * h);
    double DW = 0;
    if (r / h <= 2)
    {
        DW = alpha / h * (2 * (1 - 0.5 * r / h) * (1 - 0.5 * r / h) * (1 - 0.5 * r / h) * (1 - 0.5 * r / h) - (2 * r / h - 1) * (1 - 0.5 * r / h) * (1 - 0.5 * r / h) * (1 - 0.5 * r / h));
    }
    return DW;
}

double f_quinitc_spline(double r, double h)
{
    double alpha = 3 / (359 * M_PI * h * h * h);
    double W = 0;
    if (r / h < 1)
    {
        W = alpha * ((3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) - 6 * (2 - r / h) * (2 - r / h) * (2 - r / h) * (2 - r / h) * (2 - r / h) + 15 * (1 - r / h) * (1 - r / h) * (1 - r / h) * (1 - r / h) * (1 - r / h));
    }
    if (1 <= r / h && r / h < 2)
    {
        W = alpha * ((3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) - 6 * (2 - r / h) * (2 - r / h) * (2 - r / h) * (2 - r / h) * (2 - r / h));
    }
    if (2 <= r / h && r / h < 3)
    {
        W = alpha * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h);
    }
    return W;
}

double derive_quintic_spline(double r, double h)
{
    double alpha = 3 / (359 * M_PI * h * h * h);
    double DW = 0;
    if (r / h < 1)
    {
        DW = alpha / h * (-5 * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) + 30 * (2 - r / h) * (2 - r / h) * (2 - r / h) * (2 - r / h) - 75 * (1 - r / h) * (1 - r / h) * (1 - r / h) * (1 - r / h));
    }
    if (1 <= r / h && r / h < 2)
    {
        DW = alpha / h * (-5 * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h) + 30 * (2 - r / h) * (2 - r / h) * (2 - r / h) * (2 - r / h));
    }
    if (2 <= r / h && r / h < 3)
    {
        DW = -5 * alpha / h * (3 - r / h) * (3 - r / h) * (3 - r / h) * (3 - r / h);
    }
    return DW;
}
