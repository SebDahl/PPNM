// minimizer.h
#ifndef MINIMIZER_H
#define MINIMIZER_H

#include "matrix.h"
#include <functional>

namespace min {
pp::vector gradient(std::function<double(pp::vector)> f, pp::vector x, int mode = 0);
pp::matrix hessian(std::function<double(pp::vector)> f, pp::vector x, int mode = 0);
std::pair<pp::vector, int> newton(std::function<double(pp::vector)> f, pp::vector x, double acc = 1e-6, int mode = 0);
} // namespace min

#endif // MINIMIZER_H
