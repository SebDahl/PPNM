// minimizer.cpp
#include "minimizer.h"
#include "QR.h"
#include <cmath>
#include <stdexcept>


namespace min {
pp::vector gradient(std::function<double(pp::vector)> f, pp::vector x, int mode) {
    int n = x.size();
    pp::vector grad(n);
    double fx = f(x);
    for (int i = 0; i < n; ++i) {
        double dx = std::abs(x[i]) * std::pow(2, -26);
        x[i] += dx;
        double fx_plus = f(x);
        x[i] -= 2 * dx;
        double fx_minus = f(x);
        x[i] += dx;
        if (mode == 0)
            grad[i] = (fx_plus - fx) / dx;
        else if (mode == 1)
            grad[i] = (fx_plus - fx_minus) / (2 * dx);
        else
            throw std::runtime_error("gradient: invalid mode");
    }
    return grad;
}

pp::matrix hessian(std::function<double(pp::vector)> f, pp::vector x, int mode) {
    int n = x.size();
    pp::matrix H(n, n);
    for (int j = 0; j < n; ++j) {
        double dx = std::abs(x[j]) * std::pow(2, -26);
        x[j] += dx;
        pp::vector g_plus = gradient(f, x, mode);
        x[j] -= 2 * dx;
        pp::vector g_minus = gradient(f, x, mode);
        x[j] += dx;
        for (int i = 0; i < n; ++i) {
            if (mode == 0)
                H(i, j) = (g_plus[i] - gradient(f, x, mode)[i]) / dx;
            else if (mode == 1)
                H(i, j) = (g_plus[i] - g_minus[i]) / (2 * dx);
            else
                throw std::runtime_error("hessian: invalid mode");
        }
    }
    return H;
}

std::pair<pp::vector, int> newton(std::function<double(pp::vector)> f, pp::vector x, double acc, int mode) {
    int steps = 0;
    while (true) {
        pp::vector g = gradient(f, x, mode);
        if (g.norm() < acc) break;
        pp::matrix H = hessian(f, x, mode);
        pp::QR qr(H);
        pp::vector dx = qr.solve(-1 * g);
        double lambda = 1.0;
        while (lambda > 1.0 / 128.0 && f(x + lambda * dx) >= f(x))
            lambda /= 2;
        x += lambda * dx;
        steps++;
    }
    return {x, steps};
}
} // namespace min
