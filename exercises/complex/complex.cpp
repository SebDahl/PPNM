#include <iostream>
#include <complex>
#include <cmath>

bool approx_equal(std::complex<double> a, std::complex<double> b, double tol = 1.0e-9) {
    return std::abs(a - b) < tol;
}
void return_approx_equal(std::complex<double> a, std::complex<double> b, double tol = 1.0e-9) {
    if (approx_equal(a, b, tol)) {
        std::cout << "a and b are approximately equal." << std::endl;
    } else {
        std::cout << "a and b are not approximately equal." << std::endl;
    }
}

void print_complex(const std::complex<double>& a) {
    std::cout << "Real part: " << a.real() << "\nImaginary part: " << a.imag() << std::endl;
}

int main() {
    std::complex<double> a = std::sqrt(std::complex<double>(-1.0, 0));
    std::cout << "a = " << a << std::endl;
    print_complex(a);

    std::complex<double> b = std::log(std::complex<double>(0, 1.0));
    print_complex(b);

    std::complex<double> b_2 = std::complex<double>(0, M_PI_2);
    return_approx_equal(b, b_2);



    return 0;
}
