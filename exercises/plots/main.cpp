#include <iostream>
#include <cmath>
#include "sfuns.h"
#include <fstream>

int main() {
    // File names
    const char* erf_file = "erf_data.txt";
    const char* gamma_file = "gamma_data.txt";
    const char* lngamma_file = "lngamma_data.txt";

    // Open files for writing
    std::ofstream erf_out(erf_file);
    std::ofstream gamma_out(gamma_file);
    std::ofstream lngamma_out(lngamma_file);

    if (!erf_out || !gamma_out || !lngamma_out) {
        std::cerr << "Error opening output files!\n";
        return 1;
    }

    double x_min = -2, x_max = 2, step = 0.1;

    for (double x = x_min; x <= x_max; x += step) {
        erf_out << x << " " << sfuns::erf(x) << "\n";
        if (x > 0) // Gamma function is only valid for positive x
            gamma_out << x << " " << sfuns::fgamma(x) << "\n", lngamma_out << x << " " << sfuns::lngamma(x) << "\n";
        
    }

    erf_out.close();
    gamma_out.close();

    return 0;
}
