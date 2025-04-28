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
double corput(int N, int base){
    double q = 0, bk = 1.0/base;
    while (N > 0) {
        q += (N % base) * bk;
        N /= base;
        bk /= base;
    }
    return q;
}

std::vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
    67, 71, 73, 79, 83, 89, 97};

static pp::vector quasimc(std::function<double(pp::vector)> f, pp::vector a, pp::vector b, int N){
    int dim = a.size();
    double V=1;
    for(int i=0;i<dim;i++){V*=(b[i]-a[i]);}
    double sum=0, sum2=0;
    auto x = pp::vector(dim);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            x[j] = a[j] + (b[j] - a[j]) * (corput(i+1, primes[j]));
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



double test_function_error(pp::vector x){
    return x[0]*x[1]*x[2]/10.0;
}

double unit_circle(pp::vector x){
    return 1/(1-std::cos(x[0])*std::cos(x[1]))*1/(M_PI*M_PI);
}

double unit_circle_area = M_PI;



int main(){

    pp::vector result_list(0);
    pp::vector deviation_list(0);
    std::vector N_s = {10, 100, 1000, 10000, 100000, 1000000};

    pp::vector a = {0, 5, 10};
    pp::vector b = {5, 10, 15};
    for (double N : N_s){

        pp::vector result = plainmc(test_function_error, a, b, N);
        result_list.append(result[0]);
        deviation_list.append(result[1]);
    }
    pp::vector::write(result_list, "result_list.txt");
    pp::vector::write(deviation_list, "deviation_list.txt");
    // write N_s to file:
    std::ofstream N_s_file("N_s.txt");
    for (double N : N_s){
        N_s_file << N << std::endl;
    }
    N_s_file.close();

    


    int N = 1000000;
    a = {0, 0, 0};
    b = {M_PI, M_PI, M_PI};
    pp::vector result = plainmc(test_function, a, b, N);
    std::cout << "Integral of 1/(1-cos(x)*cos(y)*cos(z)) over [0,pi]^3" << std::endl;
    std::cout << "Result of integration: " << result[0] << std::endl;
    std::cout << "Standard deviation: " << result[1] << std::endl;
    std::cout << "Number of evaluations: " << N << std::endl;
    
    std::cout << "Test of quasimc" << std::endl;
    result = quasimc(test_function, a, b, N);
    std::cout << "Result of integration: " << result[0] << std::endl;
    std::cout << "Standard deviation: " << result[1] << std::endl;
    std::cout << "Number of evaluations: " << N << std::endl;



    return 0;
}








