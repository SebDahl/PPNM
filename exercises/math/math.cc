#include <iostream>
#include <cmath>
#include "sfuns.h"
int main() {
    double sqrt2 = std::sqrt(2);
    double fifthroot2 = std::pow(2, 0.2);
    double e_pi = std::exp(M_PI);
    double pi_e = std::pow(M_PI, M_E);
    std::cout << "sqrt(2) = " << sqrt2 << std::endl;
    std::cout << "2^(1/5) = " << fifthroot2 << std::endl;
    std::cout << "e^pi = " << e_pi << std::endl;
    std::cout << "pi^e = " << pi_e << std::endl;
    for (int i =1; i<=10; i++){
        std::cout << "gamma(" << i << ") = " << sfuns::fgamma(i) << std::endl;
    }
    return 0;
}


