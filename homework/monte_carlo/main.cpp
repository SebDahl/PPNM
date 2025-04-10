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


static pp::vector plainmc(std::function<double(pp::vector)> f, pp::vector a, pp::vector b, int N){
    int dim = a.size();
    double V=1;
    for(int i=0;i<dim;i++){V*=(b[i]-a[i]);}
    double sum=0, sum2=0;
    auto x = pp::vector(dim);
    std::random_device rd;
    std::mt19937 gen(rd());
    for (int i=0; i<N; i++){
        for (int j=0; j<dim; j++){
            x[j] = a[j] + (b[j]-a[j])*std::generate_canonical<double, 10>(gen);
        }
    double fx = f(x);
    sum += fx;
    sum2 += fx*fx;
    }
    double mean = sum/N, sigma = std::sqrt((sum2/N - mean*mean));
    pp::vector result = {mean*V, sigma*V/sqrt(N)};
    return result;
}
double test_function(pp::vector x){
    return 1/(1-std::cos(x[0])*std::cos(x[1])*std::cos(x[2]))*1/(M_PI*M_PI*M_PI);
}

int main(){
    int N = 1000000;
    pp::vector a = {0, 0, 0};
    pp::vector b = {M_PI, M_PI, M_PI};
    pp::vector result = plainmc(test_function, a, b, N);
    std::cout << "Integral of 1/(1-cos(x)*cos(y)*cos(z)) over [0,pi]^3" << std::endl;
    std::cout << "Result of integration: " << result[0] << std::endl;
    std::cout << "Standard deviation: " << result[1] << std::endl;
    std::cout << "Number of evaluations: " << N << std::endl;



    return 0;
}








