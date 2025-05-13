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

pp::vector grad_himmelblau(pp::vector x){
    double dx = 4*x[0]*(x[0]*x[0] + x[1] - 11) + 2*(x[0] + x[1]*x[1] - 7);
    double dy = 2*(x[0]*x[0] + x[1] - 11) + 4*x[1]*(x[0] + x[1]*x[1] - 7);
    return {dx, dy};
}

int main(){
    std::cout << "Check of Booth function:" << std::endl;
    std::cout << "Initial guess x = {-3, 1}" << std::endl;
    pp::vector x = {-3, 1};
    pp::vector deltax = {0.01, 0.01};
    pp::vector result = newtons_method_back(grad_booth_function, x, 1e-3, 1e-4);
    std::cout << "Result: " << result[0] <<", "<< result[1] << std::endl;
    double f_result = booth_function(result);
    std::cout << "Function value at result: " << f_result << std::endl;
    return 0;
}



