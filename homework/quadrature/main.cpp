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

int eval_count = 0;

double integrate(
    std::function<double(double)> f, /* the f from dy/dx=f(x,y) */
    double a,                        /* the lower limit of integration */
    double b,                        /* the upper limit of integration */
    double acc=0.0001,           /* absolute accuracy goal */
    double eps=0.0001,           /* relative accuracy goal */
    double f2=NAN,
    double f3=NAN // NaN indicates first call
)
{
        double h=b-a;
        
        if (std::isnan(f2)){
            eval_count -=eval_count;
            f2=f(a+2*h/6);
            eval_count += 1;
            f3=f(a+4*h/6);
            eval_count += 1;
            }
        double f1 = f(a+h/6), f4=f(a+5*h/6);
        eval_count += 2;
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

double variable_transform_integrate(
    std::function<double(double)> f,
    double a,
    double b,
    double acc=0.0001,
    double eps=0.0001,
    double f2=NAN,
    double f3=NAN // NaN indicates first call
)
{
    auto g = [=](double x){
        return f(
            (a+b)/2+(b-a)/2*std::cos(x)
        )*std::sin(x)*(b-a)/2;
    };
    return integrate(g, 0, M_PI, acc, eps);
    
}

double inf_to_inf_integrate(
    std::function<double(double)> f,
    double acc=0.0001,
    double eps=0.0001,
    double f2=NAN,
    double f3=NAN // NaN indicates first call
){
    auto g = [=](double x){
        return f(x/(1-x*x))
        *(1+x*x)/((1-x*x)*(1-x*x));
    };
    return integrate(g,-1,1, acc, eps);
}
double a_to_inf_integrate(
    std::function<double(double)> f,
    double a,
    double acc=0.0001,
    double eps=0.0001,
    double f2=NAN,
    double f3=NAN // NaN indicates first call
){
    auto g = [=](double x){
        return f(a+(1-x)/x)/(x*x);
    };
    return integrate(g,0,1,acc,eps);
}

double inf_to_b_integrate(
    std::function<double(double)> f,
    double b,
    double acc=0.0001,
    double eps=0.0001,
    double f2=NAN,
    double f3=NAN // NaN indicates first call
){
    auto g = [=](double x){
        return f(b-(1-x)/x)/(x*x);
    };
    return integrate(g,0,1,acc,eps);
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
double inverse_sqrt(double x){
    return 1/std::sqrt(x);
}

double ln_inverse_sqrt(double x){
    return std::log(x)/std::sqrt(x);
}

int main(){
    // Example usage
    double result = integrate(f, 0, M_PI);
    std::cout << "Result of integration for sin(x) [0,pi]: " << result << std::endl;
    double result2 = integrate(f2, 0, 1);   
    std::cout << "Result of integration for sqrt(x) [0,1]: " << result2 << std::endl;

    std::cout << "z = " << 0.5 << std::endl;
    std::cout << erf(0.5, 0.00001, 0.00001) << std::endl;

    pp::vector zlist = pp::vector(0);
    pp::vector erf_list = pp::vector(0);

    for (double z = -3; z <= 3; z += 0.01) {
        zlist.append(z);
        erf_list.append(erf(z, 0.001, 0.001));
    }

    pp::vector::write(erf_list, "erf_list.txt");
    pp::vector::write(zlist, "z_list.txt");

    pp::vector accs = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.00000001};
    pp::vector acc_erf = pp::vector(0);
    for (int i = 0; i < accs.size(); i++){
        double acc = accs[i];
        double result = erf(1.0, acc);
        acc_erf.append(result);
    }
    pp::vector::write(acc_erf, "acc_erf.txt");
    pp::vector::write(accs, "accs.txt");
    

    std::cout << "Test of Variable transformation quadratures" << std::endl;

    double inverse_sqrt_result = variable_transform_integrate(inverse_sqrt, 0, 1, 0.001, 0.001);
    std::cout << "Using variable transformation:" << std::endl;
    std::cout << "Result of inverse square root integration from 0 to 1 is: " << inverse_sqrt_result <<std::endl;
    std::cout << "With number of function evaluations: " << eval_count << std::endl;

    double inverse_sqrt_result_ordinary = integrate(inverse_sqrt, 0, 1, 0.001, 0.001);

    std::cout << "Using ordinary integration:" << std::endl;
    std::cout << "Result of inverse square root integration from 0 to 1 is: " << inverse_sqrt_result_ordinary <<std::endl;
    std::cout << "With number of function evaluations: " << eval_count << std::endl;


    std::cout << std::endl;
    double ln_inverse_sqrt_result = variable_transform_integrate(ln_inverse_sqrt, 0, 1, 0.001, 0.001);
    std::cout << "Result of ln inverse square root integration from 0 to 1 is: " << ln_inverse_sqrt_result <<std::endl;
    std::cout << "With number of function evaluations: " << eval_count << std::endl;

    double ln_inverse_sqrt_result_ordinary = integrate(ln_inverse_sqrt, 0, 1, 0.001, 0.001);

    std::cout << "Using ordinary integration:" << std::endl;
    std::cout << "Result of inverse square root integration from 0 to 1 is: " << ln_inverse_sqrt_result_ordinary <<std::endl;
    std::cout << "With number of function evaluations: " << eval_count << std::endl;


    double inf_gauss_result = inf_to_inf_integrate(f_gaussian, 0.001, 0.001);
    std::cout << "Integrating Gaussian from - inf to + inf: " << std::endl;
    std::cout << inf_gauss_result << std::endl;
    std::cout << "With number of evalutions: " << eval_count << std::endl;

    double a_to_inf_gauss_result = a_to_inf_integrate(f_gaussian,0, 0.001, 0.001);
    std::cout << "Integrating Gaussian from 0 to + inf: " << std::endl;
    std::cout << a_to_inf_gauss_result << std::endl;
    std::cout << "With number of evalutions: " << eval_count << std::endl;



    return 0;
}








