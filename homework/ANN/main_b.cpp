#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include "minimizer.h"

// Helper: Generate a linspace array
std::vector<double> linspace(double start, double end, int N) {
    std::vector<double> xs(N);
    double step = (end - start) / (N - 1);
    for (int i = 0; i < N; ++i)
        xs[i] = start + i * step;
    return xs;
}

// ANN with Gaussian activation
struct ann {
    int n;
    std::function<double(double)> f;
    pp::vector params, w, a, b;

    ann(int n)
        : n(n), f([](double x) { return std::exp(-x * x); }),
          params(std::vector<double>(3 * n, 1.0)) {
        update_weights_from_params();
    }

    void update_weights_from_params() {
        w = pp::vector(n);
        a = pp::vector(n);
        b = pp::vector(n);
        for (int i = 0; i < n; ++i) {
            w[i] = params[i];
            a[i] = params[i + n];
            b[i] = params[i + 2 * n];
        }
    }

    double response(double x, const pp::vector& p) const {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            double u = (x - p[i + n]) / p[i + 2 * n];
            sum += p[i] * f(u);
        }
        return sum;
    }

    double cost(const std::vector<double>& xs, double c, double y_c, double y_c_prime,
                const std::function<double(std::function<double(double)>, const std::vector<double>&)>& phi,
                const pp::vector& p,
                double dc = 1e-2) const {

        auto r = [&](double x) { return response(x, p); };
        double cost = phi(r, xs);

        double diff = r(c) - y_c;
        double rprime = (r(c + dc) - r(c - dc)) / (2 * dc);
        double diff_prime = rprime - y_c_prime;

        return cost + 100 * diff * diff + 100 * diff_prime * diff_prime;
    }

    void train(const std::vector<double>& xs, double c, double y_c, double y_c_prime,
               const std::function<double(std::function<double(double)>, const std::vector<double>&)>& phi,
               const pp::vector& guess, double acc = 1e-3) {
        auto wrapped_cost = [&](pp::vector p) {
            return cost(xs, c, y_c, y_c_prime, phi, p);
        };
        auto [opt, steps] = min::newton(wrapped_cost, guess, acc);
        params = opt;
        update_weights_from_params();
        std::cout << "Training finished in " << steps << " steps\n";
    }
};

// Oscillator parameters
constexpr double omega = 5.0;
constexpr double eta = 1.3;

// Function to compute residual integral of the damped oscillator DE
double damped_oscillator_residual(std::function<double(double)> y, const std::vector<double>& xs) {
    double total = 0;
    for (size_t i = 1; i < xs.size() - 1; ++i) {
        double h1 = xs[i] - xs[i - 1], h2 = xs[i + 1] - xs[i];
        double y_prev = y(xs[i - 1]), y_curr = y(xs[i]), y_next = y(xs[i + 1]);

        double y_prime = (y_next - y_prev) / (h1 + h2);
        double y_2prime = (y_next - 2 * y_curr + y_prev) / (0.5 * (h1 + h2) * 0.5 * (h1 + h2));

        double phi = y_2prime + eta * y_prime + omega * omega * y_curr;
        double dx = (xs[i + 1] - xs[i - 1]) / 2.0;
        total += phi * phi * dx;
    }
    return total;
}

// Analytic solution for comparison
double analytic_solution(double x) {
    return std::exp(-0.5 * eta * x) * std::cos(std::sqrt(omega * omega - 0.25 * eta * eta) * x);
}

double sign_alt(int i) { return (i % 2 == 0) ? 1.0 : -1.0; }

int main() {
    std::vector<double> xs = linspace(0, 3.0, 100);
    double c = 0, y_c = 1, y_c_prime = 0;
    int n = 8;

    ann net(n);
    pp::vector guess;
    for (int i = 0; i < n; ++i) guess.append(sign_alt(i) * 0.8 * (1 - 1.2 * double(i) / n));
    for (int i = 0; i < n; ++i) guess.append(3.0 * i / (n - 1));
    for (int i = 0; i < n; ++i) guess.append(0.2);

    std::cout << "Initial cost: " << net.cost(xs, c, y_c, y_c_prime, damped_oscillator_residual, guess) << "\n";
    net.train(xs, c, y_c, y_c_prime, damped_oscillator_residual, guess, 5);
    std::cout << "Final cost: " << net.cost(xs, c, y_c, y_c_prime, damped_oscillator_residual, net.params) << "\n";

    std::ofstream out("diffeq_values.txt");
    for (double x : xs)
        out << x << " " << net.response(x, net.params) << " " << analytic_solution(x) << "\n";
    out.close();

    return 0;
}
