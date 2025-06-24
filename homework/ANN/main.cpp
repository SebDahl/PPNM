// ann_trainer.cpp
#include <iostream>
#include <fstream>
#include <cmath>
#include "minimizer.h"

std::vector<double> linspace(double start, double end, int num_points) {
    std::vector<double> result(num_points);
    double step = (end - start) / (num_points - 1);
    for (int i = 0; i < num_points; ++i)
        result[i] = start + i * step;
    return result;
}

struct ann {
    int n;
    std::function<double(double)> f;
    pp::vector params, w, a, b;

    ann(int n) 
        : n(n), 
          f([](double x) { return x * std::exp(-x * x); }), 
          params(pp::vector(std::vector<double>(3 * n, 1.0))) 
    {
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
            double ui = (x - p[i + n]) / p[i + 2 * n];
            sum += p[i] * f(ui);
        }
        return sum;
    }

    double cost(const std::vector<double>& xs, const std::vector<double>& ys, const pp::vector& p) const {
        double sum = 0;
        for (size_t i = 0; i < xs.size(); ++i) {
            double diff = response(xs[i], p) - ys[i];
            sum += diff * diff;
        }
        return sum;
    }

    void train(const std::vector<double>& xs, const std::vector<double>& ys, const pp::vector& guess, double acc = 1e-3) {
        auto costwrap = [&](pp::vector p) { return cost(xs, ys, p); };
        auto [opt, steps] = min::newton(costwrap, guess, acc);
        params = opt;
        std::cout << "Trained with " << steps << " steps. Final cost: " << cost(xs, ys, opt) << "\n";
        update_weights_from_params();
    }

    double response_int(double x, const pp::vector& p) const {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            double ui = (x - p[i + n]) / p[i + 2 * n];
            sum += p[i] * p[i + 2 * n] * std::exp(-ui * ui);
        }
        return -0.5 * sum;
    }

    double response_deriv(double x, const pp::vector& p) const {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            double ui = (x - p[i + n]) / p[i + 2 * n];
            sum += (p[i] / p[i + 2 * n]) * (1 - 2 * ui * ui) * std::exp(-ui * ui);
        }
        return sum;
    }

    double response_2deriv(double x, const pp::vector& p) const {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            double ui = (x - p[i + n]) / p[i + 2 * n];
            sum += (p[i] / (p[i + 2 * n] * p[i + 2 * n])) * (4 * ui * ui * ui - 6 * ui) * std::exp(-ui * ui);
        }
        return sum;
    }
};

int main() {
    auto g = [](double x) { return std::cos(5 * x - 1) * std::exp(-x * x); };

    std::vector<double> xs = linspace(-1.5, 1.5, 200);
    std::vector<double> ys(xs.size());
    for (size_t i = 0; i < xs.size(); ++i)
        ys[i] = g(xs[i]);

    ann net(3);
    pp::vector guess = {-1, 1.5, -0.7, -0.7, -0.2, 0.5, 0.5, 0.5, 0.5};
    std::cout << "Initial cost: " << net.cost(xs, ys, guess) << "\n";
    net.train(xs, ys, guess, 1);

    std::ofstream valfile("values.txt");
    for (size_t i = 0; i < xs.size(); ++i) {
        valfile << xs[i] << " " << ys[i] << " "
                << net.response(xs[i], net.params) << " "
                << net.response_int(xs[i], net.params) << " "
                << net.response_deriv(xs[i], net.params) << " "
                << net.response_2deriv(xs[i], net.params) << "\n";
    }
    valfile.close();
    return 0;
}
