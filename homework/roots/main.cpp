#include <iostream>
#include <cmath>
#include "matrix.h"
#include "sfuns.h"
#include "EVD.h"
#include "QR.h"
#include "ode.h"
#include <random>
#include <fstream>
#include <functional>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>




static pp::matrix jacobian(
    std::function<pp::vector(pp::vector)> f,
    pp::vector x,
    pp::vector fx = pp::vector{},
    pp::vector deltax = pp::vector{})
{
    const int n = x.size();
    if (deltax.size() == 0) {
        deltax.resize(n);
        for (int i = 0; i < n; ++i) {
            double xi = std::abs(x[i]);
            deltax[i] = std::max(xi, 1.0) * std::pow(2.0, -26); // prevent too-small step
        }
    }

    if (fx.size() != f(x).size()) {
        fx = f(x);
    }

    const int m = fx.size();
    pp::matrix J(m, n); // support general R^n → R^m

    for (int j = 0; j < n; ++j) {
        pp::vector xj = x;
        xj[j] += deltax[j];
        pp::vector df = f(xj) - fx;
        for (int i = 0; i < m; ++i) {
            J[i, j] = df[i] / deltax[j];
        }
    }
    return J;
}


static pp::vector newtons_method_back(
    std::function<pp::vector(pp::vector)> f,
    pp::vector start, 
    double accuracy = 1e-2,
    double lambda_min = 1e-3,
    pp::vector deltax = pp::vector{0}){
        pp::vector x = start;
        pp::vector fx = f(x), z, fz;
        do{/*Newtons iterations*/
            if(fx.norm() < accuracy) break;
            pp::matrix J = jacobian(f, x, fx);
            pp::matrix::write(J, "jacobian.txt");
            auto QRJ = pp::QR(J);
            pp::vector Dx = QRJ.solve(-1*fx);
            double lambda = 1;
            do{/*line search*/
                z = x + lambda*Dx;
                fz = f(z);
                if(fz.norm() < (1 - lambda/2)*fx.norm() ) break;
                if(lambda < lambda_min) break;
                lambda *= 0.5;
            }while(true);
        x=z; fx=fz;
        }while(true);
    return x;

}

double himmelblau(pp::vector x){
    double a = x[0]*x[0] + x[1] - 11;
    double b = x[0] + x[1]*x[1] - 7;
    return a*a + b*b;
}
pp::vector grad_himmelblau(pp::vector x){
    double dx = 4*x[0]*(x[0]*x[0] + x[1] - 11) + 2*(x[0] + x[1]*x[1] - 7);
    double dy = 2*(x[0]*x[0] + x[1] - 11) + 4*x[1]*(x[0] + x[1]*x[1] - 7);
    return {dx, dy};
}

double booth_function(pp::vector x){
    double a = x[0] + 2*x[1] - 7;
    double b = 2*x[0] + x[1] - 5;
    return a*a + b*b;
}
pp::vector grad_booth_function(pp::vector x){
    double dx = 2*(x[0] + 2*x[1] - 7) + 4*(2*x[0] + x[1] - 5);
    double dy = 4*(x[0] + 2*x[1] - 7) + 2*(2*x[0] + x[1] - 5);
    return {dx, dy};
}



double rosenbrock(pp::vector x){
    double a = 1 - x[0];
    double b = x[1] - x[0]*x[0];
    return a*a + 100*b*b;
}
pp::vector grad_rosenbrock(pp::vector x){
    double dx = -2*(1-x[0]) - 400*x[0]*(x[1]-x[0]*x[0]);
    double dy = 200*(x[1]-x[0]*x[0]);
    return {dx, dy};
}

// Hydrogen atom ground state (reduced radial function)
double psi (double r) {return r*std::exp(-r);};


double rmax = 8.0, rmin = 0.1, rstep = 0.1, acc=1e-3, eps=1e-3; // All parameters
pp::vector Mfunc (pp::vector E) {
    auto schrodeq = [E](double r, pp::vector v) -> pp::vector {
        double f = v[0];
        double fprime = v[1];
        double fprimeprime = -2 * f * (E[0] + 1 / r);
        return pp::vector({fprime, fprimeprime});
    };

    std::pair<double, double> rinterval = {rmin, rmax}; // Interval for r-values
    pp::vector init({rmin-rmin*rmin, 1-2*rmin}); // Initial conditions, schrodeq(rmin)
    auto [rs, fs] = ode::driver(schrodeq, rinterval, init, rstep, acc, eps); //

    return pp::vector({fs.back()[0]}); // Return the last value of f(r), i.e. f(rmax)
};


struct TestFunction {
    std::string name;
    std::function<double(pp::vector)> f;
    std::function<pp::vector(pp::vector)> grad;
    pp::vector initial_guess;
};
int main(){
    std::vector<TestFunction> functions = {
        {"Booth", booth_function, grad_booth_function, {-3, 1}},
        {"Himmelblau", himmelblau, grad_himmelblau, {3.0, 3.0}},
        {"Rosenbrock", rosenbrock, grad_rosenbrock, {-1.2, 1.0}}
    };

    for (std::size_t i = 0; i < functions.size(); ++i) {
        const TestFunction& tf = functions[i];
        std::cout << "Check of " << tf.name << " function:" << std::endl;
        std::cout << "Initial guess x = {" << tf.initial_guess[0] << ", " << tf.initial_guess[1] << "}" << std::endl;
        pp::vector result = newtons_method_back(tf.grad, tf.initial_guess, 1e-3, 1e-4);
        std::cout << "Result: " << result[0] << ", " << result[1] << std::endl;
        double f_result = tf.f(result);
        std::cout << "Function value at result: " << f_result << std::endl;
        std::cout << std::endl;
    }

    // Test the Mfunc function
    pp::vector E = {0.5};
    pp::vector result = Mfunc(E);
    std::cout << "Result of Mfunc: " << result[0] << std::endl;
    std::cout << "Expected value: " << std::exp(-2.0) << std::endl;
    std::cout << "Difference: " << std::abs(result[0] - std::exp(-2.0)) << std::endl;

    return 0;
}



