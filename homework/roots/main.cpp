// Newton's method with numerical Jacobian and hydrogen atom root finder
#include <iostream>
#include <cmath>
#include <functional>
#include <string>
#include <vector>
#include <fstream>
#include "matrix.h"
#include "QR.h"
#include "ode.h"

// Numerical Jacobian
static pp::matrix jacobian(
    std::function<pp::vector(pp::vector)> f,
    pp::vector x,
    pp::vector fx = pp::vector{},
    pp::vector deltax = pp::vector{})
{
    int n = x.size();
    if (deltax.size() != n) {
        deltax.resize(n);
        for (int i = 0; i < n; ++i)
            deltax[i] = std::max(std::abs(x[i]), 1.0) * std::pow(2.0, -26);
    }
    if (fx.size() != f(x).size()) fx = f(x);
    int m = fx.size();
    pp::matrix J(m, n);
    for (int j = 0; j < n; ++j) {
        x[j] += deltax[j];
        pp::vector df = f(x) - fx;
        for (int i = 0; i < m; ++i)
            J[i, j] = df[i] / deltax[j];
        x[j] -= deltax[j];
    }
    return J;
}

// Newton's method with backtracking line search
static pp::vector newton(
    std::function<pp::vector(pp::vector)> f,
    pp::vector start,
    double acc = 1e-2,
    double lambda_min = 1e-3,
    pp::vector deltax = pp::vector{})
{
    pp::vector x = start.copy(), fx = f(x), z, fz;
    while (true) {
        if (fx.norm() < acc) break;
        pp::matrix J = jacobian(f, x, fx, deltax);
        pp::QR qrJ(J);
        pp::vector Dx = qrJ.solve(-1 * fx);
        double lambda = 1.0;
        do {
            z = x + lambda * Dx;
            fz = f(z);
            if (fz.norm() < (1 - lambda / 2) * fx.norm()) break;
            if (lambda < lambda_min) break;
            lambda *= 0.5;
        } while (true);
        x = z; fx = fz;
    }
    return x;
}

// Hydrogen Schrödinger M-function
pp::vector Mfunc(pp::vector E, double rmin, double rmax, double h, double acc, double eps, std::vector<double>& rs_out, std::vector<double>& fs_out) {
    auto schrodeq = [E](double r, const pp::vector& v) -> pp::vector {
        double f = v[0];
        double fprime = v[1];
        double fprimeprime = -2 * f * (E[0] + 1.0 / r);
        return {fprime, fprimeprime};
    };
    pp::vector init = {rmin - rmin * rmin, 1 - 2 * rmin};
    auto [rvals, fvals] = ode::driver(schrodeq, {rmin, rmax}, init, h, acc, eps);
    rs_out.clear(); fs_out.clear();
    for (int i = 0; i < rvals.size(); ++i) {
        rs_out.push_back(rvals[i]);
        fs_out.push_back(fvals[i]);
    }
    return {fvals[fvals.size() - 1]};
}

int main() {
    // Hydrogen shooting method root with convergence analysis
    std::vector<std::pair<std::string, std::tuple<double, double, double, double>>> tests = {
        {"baseline", {0.1, 8.0, 0.1, 1e-6}},
        {"rmax_10", {0.1, 10.0, 0.1, 1e-6}},
        {"rmin_0.01", {0.01, 8.0, 0.1, 1e-6}},
        {"acc_1e-8", {0.1, 8.0, 0.1, 1e-8}},
        {"eps_1e-8", {0.1, 8.0, 0.1, 1e-6}},
    };

    for (auto& [label, params] : tests) {
        double rmin, rmax, h, acc;
        std::tie(rmin, rmax, h, acc) = params;
        double eps = (label == "eps_1e-8" ? 1e-8 : 1e-6);
        std::vector<double> rs, fs;
        auto M = [&](pp::vector E) { return Mfunc(E, rmin, rmax, h, acc, eps, rs, fs); };
        pp::vector E0 = newton(M, {-0.7}, 1e-8);

        std::cout << "Test: " << label << "\n";
        std::cout << "  Found E0 = " << E0[0] << ", Error = " << std::abs(E0[0] + 0.5) << "\n";

        std::ofstream wfout("wavefunction_" + label + ".txt");
        for (size_t i = 0; i < rs.size(); ++i) {
            double exact = rs[i] * std::exp(-rs[i]);
            wfout << rs[i] << " " << fs[i] << " " << exact << "\n";
        }
        wfout.close();
    }

    return 0;
}
