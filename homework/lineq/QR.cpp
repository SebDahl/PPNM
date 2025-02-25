#include "QR.h"
#include <string>
#include <algorithm>
#include <stdexcept> // For exception handling
#include <cassert>   // For assert

QR::QR(const matrix<double>& A) {
    int m = A.sizerow();
    int n = A.sizecol();
    assert(m >= n && "Matrix must have more rows than columns (m >= n)");

    matrix<double> Q(m, n);  // Q is m × n
    matrix<double> R(n, n);  // R is n × n (upper triangular)
    matrix<double> A_copy = A;

    for (int i = 0; i < n; i++) {
        vector<double> a = A_copy.get_col(i);
        double a_norm = a.norm();
        if (std::abs(a_norm) < 1e-10) {
            throw std::runtime_error("Matrix is rank-deficient: cannot compute QR.");
        }

        R(i, i) = a_norm;
        vector<double> q = a / a_norm;
        Q.set_col(i, q);

        for (int j = i + 1; j < m; j++) {
            vector<double> b = A_copy.get_col(j);
            double dot = q.dot(b);
            R(i, j) = dot;
            vector<double> new_b = b - q * dot;
            A_copy.set_col(j, new_b);
        }
    }
}


QR::solve(const std::vector<double>& b){
    int n_b = b.size();
    assert(n_b == Q.sizerow() && "Matrix Q and vector b must have the same number of rows.");
    


}


