#include <iostream>
#include <cmath>

double f_gaussian(double r, double h) {
    double alpha = 1/(pow(M_PI,3/2)*pow(h,3));
    double W = alpha * exp(-pow(r/h, 2));
    return W;
}

double deriv_gaussian(double r, double h) {
    double alpha = 1/(pow(M_PI,3/2)*pow(h,3));
    double DW = alpha/h * (-2*r/h*exp(-pow(r/h, 2)));
    return DW;
}

double f_bell(double r, double h){
    double alpha = 105/(16*M_PI*pow(h,3));
    double W = 0;
    if (r/h <= 1){
        W = (1+3*r/h)*pow(1-r/h,3);
    }
    return W;
}

double derive_bell(double r, double h){
    double alpha = 105/(16*M_PI*pow(h,3));
    double DW = 0;
    if (r/h <= 1){
        DW = alpha/h * 3 * (pow(1-r/h,3)-(1+3*r/h)*pow(1-r/h,3));
    }
    return DW;
}

double f_cubic_spline(double r, double h)
{   double alpha = 3/(2*M_PI*pow(h,3));
    double W = 0;
    if(r/h < 1){
        W = 3/2 - pow(r/h,2) + 1/2 * pow(r/h,3);
    }
    if(1<= r/h <2){
        W = 1/6 * pow(1-r/h,3);
    }
    return W;
}

double derive_cubic_spline(double r, double h)
{   double alpha = 3/(2*M_PI*pow(h,3));
    double DW = 0;
    }
    if(1<= r/h <2){
        DW = alpha/h *(3/2*pow(r/h,2) - 2*r/h);
    if(r/h < 1){
        DW = alpha/h * (-1/2*pow(2-r/h,2));
    }
    return W;
}

double f_quadratic(double r, double h){
    double alpha = 5/(4*M_PI*pow(h,3));
    double W = 0;
    if(r/h <=2){
        W = alpha *(3/16*pow(r/h,2)-3/4*r/h + 3/4);
    }
    return W;
}

double derive_quadratic(double r, double h){
    double alpha = 5/(4*M_PI*pow(h,3));
    double DW = 0;
    if (r/h <=2){
        DW = alpha/h *(3/8*r/h-3/4);
    }
    return DW;
}

double f_quintic(double r, double h){
    double alpha = 21/(16*M_PI*pow(h,3));
    double W = 0;
    if(r/h<=2){
        W = alpha*(pow(1-1/2*r/h,4)*(2*r/h+1));
    }

    return W;
}

double derive_quintic(double r, double h){
    double alpha = 21/(16*M_PI*pow(h,3));
    double DW = 0;
    if(r/h <= 2){
        DW = alpha/h *(-5*r/h*pow(1-1/2*r/h,3));

    }
    return DW;
}

double f_quinitc_spline(double r, double h){
    double alpha = 3/(359*M_PI*pow(h,3));
    double W = 0;
    if(r/h <1){
        W = alpha*(pow(3-r/h,5)-6*pow(2-r/h,5)+15*pow(1-r/h,5));
    }
    if(1<= r/h < 2){
        W = alpha*(pow(3-r/h,5)-6*pow(2-r/h,5));
    }
    if(2<=r/h<3){
        W = alpha*pow(3-r/h,5);
    }
    return W;

}

double derive_quintic_spline(double r, double h){
    double alpha = 3/(359*M_PI*pow(h,3));
    double DW = 0;
    if(r/h <1){
        DW = alpha/h*(-5*pow(3-r/h,4)+30*pow(2-r/h,4)-75*pow(1-r/h,4));
    }
    if(1<= r/h < 2){
       DW = alpha/h*(-5*pow(3-r/h,4)+30*pow(2-r/h,4));
    }
    if(2<=r/h<3){
        DW = -5*alpha/h*pow(3-r/h,4);
    }
    return DW;   
}
