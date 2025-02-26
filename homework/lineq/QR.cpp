#include "QR.h"
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cassert>

namespace pp {

QR::QR(const matrix& A) {
    int m = A.sizerow();
    int n = A.sizecol();
    assert(m >= n && "Matrix must have more rows than columns (m >= n)");

    Q = matrix(m, n);  // Q is m × n
    R = matrix(n, n);  // R is n × n (upper triangular)
    matrix A_copy = A;

    for (int i = 0; i < n; i++) {
        vector a = A_copy[i];  
        double a_norm = 0.0;

        for (int j = 0; j < a.size(); j++) {
            a_norm += a[j] * a[j];
        }
        a_norm = std::sqrt(a_norm);

        if (std::abs(a_norm) < 1e-10) {
            throw std::runtime_error("Matrix is rank-deficient: cannot compute QR.");
        }

        R(i, i) = a_norm;
        vector q = a / a_norm;
        Q[i] = q;  // Assuming `set_col(i, q)` should be replaced with `Q[i] = q`

        for (int j = i + 1; j < n; j++) {
            vector b = A_copy[j];  // Assuming `get_col(j)` should be replaced with `A_copy[j]`
            double dot = 0.0;

            for (int k = 0; k < b.size(); k++) {
                dot += q[k] * b[k];
            }

            R(i, j) = dot;
            vector new_b = b - q * dot;
            A_copy[j] = new_b;  // Assuming `set_col(j, new_b)` should be replaced with `A_copy[j] = new_b`
        }
    }
}

vector QR::solve(const vector& b) const {
    int n = Q.sizecol();
    int m = Q.sizerow();
    assert(b.size() == m && "Matrix Q and vector b must have the same number of rows.");

    // Compute y = Q^T * b
    vector y(n);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < m; j++) {
            sum += Q[j][i] * b[j];
        }
        y[i] = sum;
    }

    // Solve Rx = y using back-substitution
    vector x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int k = i + 1; k < n; k++) {
            sum += R(i, k) * x[k];
        }
        if (std::abs(R(i, i)) < 1e-10) {
            throw std::runtime_error("Singular matrix: Cannot perform back-substitution.");
        }
        x[i] = (y[i] - sum) / R(i, i);
    }
    return x;
}

double QR::det() const {
    int n = R.sizerow();
    double determinant = 1.0;
    for (int i = 0; i < n; i++) {
        determinant *= R(i, i);
    }
    return determinant;
}

}  // namespace pp
