#include <iostream>
#include <cmath>
#include "matrix.h"
#include "EVD.h"
#include "QR.h"
#include <random>
#include <fstream>
#include <functional>
#include <string>
#include <stdexcept>
#include <vector>

int binsearch(pp::vector x, double z) {
    if (z < x[0] || z > x[x.size()-1]) throw std::runtime_error("binsearch: bad z");
    int i=0, j=x.size()-1;
    while (j-i>1) {
        int mid=(i+j)/2;
        if (z>x[mid]) i=mid; else j=mid;
    }
    return i;
}
struct qspline {
    pp::vector x, y, b, c; // spline coefficients: f_i(x) = y_i + b_i*(x - x_i) + c_i*(x - x_i)^2

    qspline(const pp::vector& xs, const pp::vector& ys) : x(xs), y(ys), b(xs.size()-1), c(xs.size()-1) {
        int n = x.size();
        pp::vector dx(n-1), dy(n-1);
        for (int i = 0; i < n-1; i++) {
            dx[i] = x[i+1] - x[i];
            dy[i] = y[i+1] - y[i];
        }
        for (int i = 0; i < n-1; i++) {
            c[i] = (dy[i+1] / dx[i+1] - dy[i] / dx[i]) / (dx[i] + dx[i+1]);
        }
        c[0] = c[n-2] = 0; // natural boundary condition
        for (int i = 0; i < n-1; i++) {
            b[i] = dy[i]/dx[i] - c[i]*dx[i];
        }
    }

    double evaluate(double z) const {
        int i = binsearch(x, z);
        double dx = z - x[i];
        return y[i] + b[i]*dx + c[i]*dx*dx;
    }

    double derivative(double z) const {
        int i = binsearch(x, z);
        double dx = z - x[i];
        return b[i] + 2*c[i]*dx;
    }

    double integral(double z) const {
        int i = binsearch(x, z);
        double integral = 0.0;
        for (int j = 0; j < i; j++) {
            double dx = x[j+1] - x[j];
            integral += y[j]*dx + b[j]*dx*dx/2 + c[j]*dx*dx*dx/3;
        }
        double dx = z - x[i];
        integral += y[i]*dx + b[i]*dx*dx/2 + c[i]*dx*dx*dx/3;
        return integral;
    }
};



double linterp(pp::vector x, pp::vector y, double z) {
    int i = binsearch(x, z);
    double dx = x[i+1] - x[i];
    if (!(dx > 0)) throw std::runtime_error("uups...");
    double dy = y[i+1] - y[i];
    return y[i] + dy/dx*(z - x[i]);
}

double linterpIntegrate(pp::vector x, pp::vector y, double z) {
    if (z < x[0] || z > x[x.size() - 1]) throw std::runtime_error("linterpIntegrate: z out of bounds");
    double integral = 0.0;
    int i;
    for(i = 0; x[i+1] < z; i++) {
        double dx = x[i+1] - x[i];
        double dy = y[i+1] - y[i];
        double slope = dy / dx;
        integral += y[i]*dx + 0.5*slope*dx*dx;
    }
    double dx = z - x[i];
    double dy = y[i+1] - y[i];
    double slope = dy / (x[i+1] - x[i]);
    integral += y[i]*dx + 0.5*slope*dx*dx;
    return integral;
}

int main() {
    int N = 10;
    pp::vector x_i(N), y_i(N);
    for (int i = 0; i < N; i++) {
        x_i[i] = i;
        y_i[i] = std::sin(x_i[i]);
    }

    pp::vector::write(x_i, "x_i.txt");
    pp::vector::write(y_i, "y_i.txt");

    pp::vector x_plot, y_linterp, y_lintegral, y_qspline, y_qderiv, y_qintegral;
    qspline q(x_i, y_i);

    for (double z = 0; z < x_i[N-1]; z += 0.1) {
        x_plot.append(z);
        y_linterp.append(linterp(x_i, y_i, z));
        y_lintegral.append(linterpIntegrate(x_i, y_i, z));
        y_qspline.append(q.evaluate(z));
        y_qderiv.append(q.derivative(z));
        y_qintegral.append(q.integral(z));
    }

    pp::vector::write(x_plot, "x_plot.txt");
    pp::vector::write(y_linterp, "y_linterp.txt");
    pp::vector::write(y_lintegral, "y_lintegral.txt");
    pp::vector::write(y_qspline, "y_qspline.txt");
    pp::vector::write(y_qderiv, "y_qderiv.txt");
    pp::vector::write(y_qintegral, "y_qintegral.txt");

    return 0;
}
