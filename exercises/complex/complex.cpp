#include <iostream>
#include <complex>
#include <cmath>

bool approx_equal(std::complex<double> a, std::complex<double> b, double tol = 1.0e-9) {
    return std::abs(a - b) < tol;
}

void return_approx_equal(std::complex<double> a, std::complex<double> b, double tol = 1.0e-9) {
    std::cout << "Expected: " << b << "\nComputed: " << a << std::endl;
    if (approx_equal(a, b, tol)) {
        std::cout << "a and b are approximately equal.\n" << std::endl;
    } else {
        std::cout << "a and b are NOT approximately equal.\n" << std::endl;
    }
}

int main() {
    using namespace std;

    // 1. sqrt(-1) = ±i
    std::cout << "Checking sqrt(-1) = ±i" << std::endl;
    complex<double> a = sqrt(complex<double>(-1.0, 0));
    return_approx_equal(a, complex<double>(0, 1));

    // 2. sqrt(i) = e^(½ln(i)) = e^(iπ/4) = cos(π/4) + i sin(π/4)
    std::cout << "Checking sqrt(i) = e^(½ln(i)) = e^(iπ/4) = cos(π/4) + i sin(π/4)" << std::endl;
    complex<double> sqrt_i = sqrt(complex<double>(0, 1));
    complex<double> expected_sqrt_i = polar(1.0 / sqrt(2.0), M_PI / 4);
    return_approx_equal(sqrt_i, expected_sqrt_i);

    // 3. e^i = cos(1) + i sin(1)
    std::cout << "Checking e^i = cos(1) + i sin(1)" << std::endl;
    complex<double> e_to_i = exp(complex<double>(0, 1));
    complex<double> expected_e_to_i = complex<double>(cos(1.0), sin(1.0));
    return_approx_equal(e_to_i, expected_e_to_i);

    // 4. e^{iπ} = -1 (Euler's identity)
    std::cout << "Checking e^{iπ} = -1 (Euler's identity)" << std::endl;
    complex<double> e_to_ipi = exp(complex<double>(0, M_PI));
    return_approx_equal(e_to_ipi, complex<double>(-1, 0));

    // 5. i^i = e^{i ln(i)} = e^{-π/2} ≈ 0.207879...
    std::cout << "Checking i^i = e^{i ln(i)} = e^{-π/2} ≈ 0.207879..." << std::endl;
    complex<double> i_pow_i = pow(complex<double>(0, 1), complex<double>(0, 1));
    double expected_real = exp(-M_PI / 2);
    return_approx_equal(i_pow_i, complex<double>(expected_real, 0));

    // 6. ln(i) = ln(e^{iπ/2}) = iπ/2
    complex<double> ln_i = log(complex<double>(0, 1));
    complex<double> expected_ln_i = complex<double>(0, M_PI / 2);
    return_approx_equal(ln_i, complex<double>(0, M_PI_2));

    // 7. sin(iπ) = i sinh(π)
    std::cout << "Checking sin(iπ) = i sinh(π)" << std::endl;
    complex<double> sin_ipi = sin(complex<double>(0, M_PI));
    complex<double> expected_sin_ipi = complex<double>(0, sinh(M_PI));
    return_approx_equal(sin_ipi, expected_sin_ipi);

    return 0;
}
