#include <iostream>
#include <cmath>
#include "matrix.h"
#include "sfuns.h"
#include "EVD.h"
#include "QR.h"
#include <random>
#include <fstream>
#include <functional>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>


double integrate(
    std::function<double(double)> f, /* the f from dy/dx=f(x,y) */
    double a,                        /* the lower limit of integration */
    double b,                        /* the upper limit of integration */
    double acc=0.001,           /* absolute accuracy goal */
    double eps=0.0001,           /* relative accuracy goal */
    double f2=NAN,
    double f3=NAN // NaN indicates first call
)
{
        double h=b-a;
        if (std::isnan(f2)){
            f2=f(a+2*h/6);
            f3=f(a+4*h/6); 
            }
        double f1 = f(a+h/6), f4=f(a+5*h/6);
        double Q = (2*f1+f2+f3+2*f4)/6*h;
        double q = (f1+f2+f3+f4)/4*h;
        double err = std::abs(Q-q);
        if (err <= acc + eps*std::abs(Q)){
            return Q;
        }
        else{
            double Q1 = integrate(f, a, (a+b)/2, acc/sqrt(2), eps, f1, f2);
            double Q2 = integrate(f, (a+b)/2, b, acc/sqrt(2), eps, f3, f4);
            return Q1+Q2;
        }
}

    double f(double x) {
        return std::sin(x); // Example function
    }

    double f2(double x) {
        return std::sqrt(x); // Example function
    }
    double f_gaussian(double x) {
        return std::exp(-x*x);
    }

double erf(double z, double acc=0.001, double eps=0.001){
    if(z<0){
        return -erf(-z, acc, eps);
    }
    if (0<=z<=1){
        return 2./std::sqrt(M_PI)*integrate(f_gaussian, 0, z, acc, eps);
    }
    if (1<z){
        auto f_t = [z](double t) {
            return std::exp(-(z + (1. - t)/t)*(z + (1 - t)/t)) / (t * t);
        };
        return 1. - 2/std::sqrt(M_PI)*integrate(f_t, 0, 1, acc, eps);
    }
    
}

int main(){
    // Example usage
    double result = integrate(f, 0, M_PI);
    std::cout << "Result of integration for sin(x) [0,pi]: " << result << std::endl;
    double result2 = integrate(f2, 0, 1);
    std::cout << "Result of integration for sqrt(x) [0,1]: " << result2 << std::endl;

    std::cout << "z = " << 0.5 << std::endl;
    std::cout << erf(0.5, 0.001, 0.001) << std::endl;

    pp::vector zlist = pp::vector(0);
    pp::vector erf_list = pp::vector(0);

    for (double z = -3; z <= 3; z += 0.01) {
        zlist.append(z);
        erf_list.append(erf(z, 0.001, 0.001));
    }

    pp::vector::write(erf_list, "erf_list.txt");
    pp::vector::write(zlist, "z_list.txt");

    pp::vector accs = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.0000000000001};
    pp::vector acc_erf = pp::vector(0);
    for (int i = 0; i < accs.size(); i++){
        double acc = accs[i];
        double result = erf(1.0, acc);
        acc_erf.append(result);
    }
    pp::vector::write(acc_erf, "acc_erf.txt");
    pp::vector::write(accs, "accs.txt");




    return 0;
}








