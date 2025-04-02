#include <iostream>
#include <cmath>
#include "matrix.h"
#include "EVD.h"
#include "QR.h"
#include <random>
#include <fstream>
#include <functional>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>


struct RKStep12Result {
    pp::vector yh;
    pp::vector delta_y;
};

RKStep12Result rkstep12(
    std::function<pp::vector(double, const pp::vector&)> f, /* the f from dy/dx=f(x,y) */
    double x,                                               /* the current value of the variable */
    pp::vector y,                                           /* the current value y(x) of the sought function */
    double h                                                /* the step to be taken */
)
{
    pp::vector k0 = f(x, y);                 /* embedded lower order formula (Euler) */
    pp::vector k1 = f(x+h/2, y+k0*(h/2));    /* higher order formula (midpoint) */
    pp::vector yh = y+k1*h;                  /* y(x+h) estimate */
    pp::vector delta_y = (k1-k0)*h;          /* error estimate */
    return {yh, delta_y};
};

struct driver_result {
    pp::vector xlist;
    pp::vector ylist;
};

driver_result driver(
    std::function<pp::vector(double, const pp::vector&)> f, /* the f from dy/dx=f(x,y) */
    std::pair<double, double> interval,                              /* the interval (a,b) */
    pp::vector yinit,                                       /* the initial value y(a) */
    double h=0.125,                                        /* initial step size */
    double acc=0.1,                                        /* absolute accuracy goal */
    double eps=0.01                                     /* relative accuracy goal */
){
    double a = interval.first; double b = interval.second; double x = a; pp::vector y = yinit.copy();
    pp::vector xlist = pp::vector(0); xlist.append(x);
    pp::vector ylist = pp::vector(0); ylist.append(y[0]);
    do{
        if(x>=b) return {xlist, ylist};
        if(x+h>b) h=b-x;
        auto [yh, delta_y] = rkstep12(f, x, y, h);
        double tol = (acc+eps*yh.norm())*std::sqrt(h/(b-a));
        double err = delta_y.norm();
        if(err<=tol){//accept step
            x+=h; y=yh;
            xlist.append(x);
            ylist.append(y[0]);
            }
        if(err>0) h*=std::min(0.95*std::pow(tol/err,0.25), 2.); //readjust step size
        else h*=2;
    }while(true);
}//driver


int main(){
// Some main function
    pp::vector yinit(2); yinit[0] = 1; yinit[1] = 3;
    auto f = [](double x, const pp::vector& y) -> pp::vector {
        pp::vector dy(1);
        double rate = 0.5;  // Growth rate (positive for growth, negative for decay)
        dy[0] = rate * y[0];
        return dy;
    };
    
    driver_result result = driver(f, {0, 10}, yinit);
    pp::vector::write(result.xlist, "xlist.txt");
    pp::vector::write(result.ylist, "ylist.txt");

    return 0;
}








