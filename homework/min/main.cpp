// Optimized and adapted version using pp::matrix and pp::vector
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include "matrix.h"
#include "QR.h"

// --- Utility ---

std::vector<double> linspace(double start, double end, int num_points) {
    std::vector<double> result(num_points);
    double step = (end - start) / (num_points - 1);
    for (int i = 0; i < num_points; ++i)
        result[i] = start + i * step;
    return result;
}

void load_higgs_data(const std::string& filename, std::vector<double>& E, std::vector<double>& sigma, std::vector<double>& delta_sigma) {
    std::ifstream infile(filename);
    double e, s, d;
    while (infile >> e >> s >> d) {
        E.push_back(e);
        sigma.push_back(s);
        delta_sigma.push_back(d);
    }
}

// --- Physics Models ---

// Rosenbrock function
double rbvalley(const pp::vector& v) {
    return std::pow(1 - v[0], 2) + 100 * std::pow(v[1] - v[0] * v[0], 2);
}

// Himmelblau function
double himblau(const pp::vector& v) {
    return std::pow(v[0]*v[0] + v[1] - 11, 2) + std::pow(v[0] + v[1]*v[1] - 7, 2);
}

// Breit-Wigner
double BW(double E, const pp::vector& p) {
    double m = p[0], g = p[1], A = p[2];
    return A / (std::pow(E - m, 2) + g * g / 4);
}

// Chi-squared
double chi2(const std::vector<double>& E, const std::vector<double>& sigma, const std::vector<double>& dsigma, const pp::vector& p) {
    double chi = 0.0;
    for (size_t i = 0; i < E.size(); ++i) {
        double diff = sigma[i] - BW(E[i], p);
        chi += (diff / dsigma[i]) * (diff / dsigma[i]);
    }
    return chi;
}

// --- Derivative tools ---

pp::vector gradient(std::function<double(pp::vector)> f, pp::vector x, int mode = 0) {
    double fx = f(x);
    int n = x.size();
    pp::vector grad(n);
    for (int i = 0; i < n; ++i) {
        double dx = std::abs(x[i]) * std::pow(2, -26);
        x[i] += dx;
        double fx_plus = f(x);
        x[i] -= 2 * dx;
        double fx_minus = f(x);
        x[i] += dx;
        grad[i] = (mode == 0) ? (fx_plus - fx) / dx : (fx_plus - fx_minus) / (2 * dx);
    }
    return grad;
}

pp::matrix hessian(std::function<double(pp::vector)> f, pp::vector x, int mode = 0) {
    int n = x.size();
    pp::matrix H(n, n);
    pp::vector g0 = gradient(f, x, mode);
    for (int j = 0; j < n; ++j) {
        double dx = std::abs(x[j]) * std::pow(2, -26);
        x[j] += dx;
        pp::vector g_plus = gradient(f, x, mode);
        x[j] -= 2 * dx;
        pp::vector g_minus = gradient(f, x, mode);
        x[j] += dx;
        for (int i = 0; i < n; ++i) {
            H(i, j) = (mode == 0) ? (g_plus[i] - g0[i]) / dx : (g_plus[i] - g_minus[i]) / (2 * dx);
        }
    }
    return H;
}

// --- Newton optimizer ---

std::pair<pp::vector, int> newton(std::function<double(pp::vector)> f, pp::vector x, double acc = 1e-6, int mode = 0) {
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

// --- Main ---

int main() {
    double acc = 1e-6;
    std::vector<int> modes = {0, 1};
    std::vector<double> E, sigma, dsigma;
    load_higgs_data("higgs_data.data", E, sigma, dsigma);

    for (int mode : modes) {
        std::cout << "--- Mode: " << (mode == 0 ? "Forward" : "Central") << " difference ---\n";

        auto [rbmin, rbsteps] = newton(rbvalley, {4, 5}, acc, mode);
        std::cout << "Rosenbrock min: "; rbmin.print();
        std::cout << "Steps: " << rbsteps << "\n\n";

        auto [hbmin, hbsteps] = newton(himblau, {2.5, 3.5}, acc, mode);
        std::cout << "Himmelblau min: "; hbmin.print();
        std::cout << "Steps: " << hbsteps << "\n\n";

        auto chi2_func = [&](pp::vector p) { return chi2(E, sigma, dsigma, p); };
        auto [higgsp, hsteps] = newton(chi2_func, {125, 5, 25}, acc, mode);

        std::cout << "Higgs BW optimal parameters: "; higgsp.print();
        std::cout << "Steps: " << hsteps << ", chi2 = " << chi2(E, sigma, dsigma, higgsp) << "\n\n";

        std::ofstream fit("higgs_fit.txt");
        for (double e : linspace(100, 160, 500))
            fit << e << " " << BW(e, higgsp) << "\n";
    }

    return 0;
}
