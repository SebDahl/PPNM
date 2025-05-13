#ifndef ODE_H
#define ODE_H

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

namespace ode {

    // Structure to store result of rkstep12
    struct RKStep12Result {
        pp::vector yh;
        pp::vector delta_y;
    };

    // One step of the second-order Runge-Kutta method (midpoint method)
    RKStep12Result rkstep12(
        std::function<pp::vector(double, const pp::vector&)> f,
        double x,
        pp::vector y,
        double h
    );

    // Structure to store result of the ODE driver
    struct driver_result {
        pp::vector xlist;
        pp::vector ylist;
    };

    // Adaptive step size ODE driver
    driver_result driver(
        std::function<pp::vector(double, const pp::vector&)> f,
        std::pair<double, double> interval,
        pp::vector yinit,
        double h = 0.125,
        double acc = 0.1,
        double eps = 0.01
    );

} // namespace ode

#endif // ODE_H
