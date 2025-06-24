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
#include <algorithm>

#include "ode.h"

namespace ode {

RKStep12Result rkstep12(
    std::function<pp::vector(double, const pp::vector&)> f,
    double x,
    pp::vector y,
    double h
) {
    pp::vector k0 = f(x, y);
    pp::vector k1 = f(x + h / 2, y + k0 * (h / 2));
    pp::vector yh = y + k1 * h;
    pp::vector delta_y = (k1 - k0) * h;
    return {yh, delta_y};
}

driver_result driver(
    std::function<pp::vector(double, const pp::vector&)> f,
    std::pair<double, double> interval,
    pp::vector yinit,
    double h,
    double acc,
    double eps
) {
    double a = interval.first, b = interval.second, x = a;
    pp::vector y = yinit.copy();
    pp::vector xlist = pp::vector(0); xlist.append(x);
    pp::vector ylist = pp::vector(0); ylist.append(y[0]);

    do {
        if (x >= b) return {xlist, ylist};
        if (x + h > b) h = b - x;

        auto [yh, delta_y] = rkstep12(f, x, y, h);
        double tol = (acc + eps * yh.norm()) * std::sqrt(h / (b - a));
        double err = delta_y.norm();

        if (err <= tol) {
            x += h; y = yh;
            xlist.append(x); ylist.append(y[0]);
        }

        h *= (err > 0) ? std::min(0.95 * std::pow(tol / err, 0.25), 2.0) : 2;
    } while (true);
}

}














