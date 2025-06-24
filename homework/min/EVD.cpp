#include "QR.h"
#include "EVD.h"
#include "matrix.h"
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cassert>

namespace pp {
    void EVD::timesJ(pp::matrix& A, int p, int q, double theta){
        double c = std::cos(theta), s = std::sin(theta);
        for (int i = 0; i < A.sizerow(); i++){
            double A_ip = A(i, p), A_iq = A(i, q);
            A(i, p) = c * A_ip - s * A_iq;
            A(i, q) = s * A_ip + c * A_iq;

        }
    }
    void EVD::Jtimes(pp::matrix& A, int p, int q, double theta){
        double c = std::cos(theta), s = std::sin(theta);
        for (int j = 0; j < A.sizecol(); j++){
            double A_pj = A(p, j), A_qj = A(q, j);
            A(p, j) = c * A_pj + s * A_qj;
            A(q, j) = -s * A_pj + c * A_qj;
        }

    }
    EVD::EVD(const pp::matrix& A) {
        int n = A.sizerow();
        assert(A.sizecol() == n && "Matrix must be square (n x n)");
    
        w = pp::vector(n);  // Eigenvalues
        V = pp::matrix::identity(n); // Start with identity for eigenvectors
        D = pp::matrix::identity(n);
    
        pp::matrix A_copy = A;  // Copy of A to be modified
        cyclic_jacobi(A_copy);
    
        // Copy diagonal elements to w
        for (int i = 0; i < n; i++)
            w[i] = A_copy(i, i);
        for (int i = 0; i < n; i++)
            D(i, i) = w[i];
    }


// In place Jacobian matrix multiplication:


void EVD::cyclic_jacobi(pp::matrix& A) {
    int n = A.sizerow();
    bool changed = true;

    while (changed) {
        changed = false;
        for (int p = 0; p < n - 1; p++) {
            for (int q = p + 1; q < n; q++) {
                double apq = A(p, q), app = A(p, p), aqq = A(q, q);
                double theta = 0.5 * std::atan2(2 * apq, aqq - app);
                double c = std::cos(theta), s = std::sin(theta);

                double new_app = c * c * app - 2 * s * c * apq + s * s * aqq;
                double new_aqq = s * s * app + 2 * s * c * apq + c * c * aqq;

                // if (std::abs(new_app - app) > 1e-13 || std::abs(new_aqq - aqq) > 1e-13) {
                if (new_app != app|| new_aqq != aqq) {

                    changed = true;
                    timesJ(A, p, q, theta);
                    Jtimes(A, p, q, -theta);
                    timesJ(V, p, q, theta); // Update eigenvectors
                }
            }
        }
    }
}
// Compute determinant (product of eigenvalues)
double EVD::det() const {
    double determinant = 1.0;
    for (int i = 0; i < w.size(); i++)
        determinant *= w[i];
    return determinant;
}
pp::matrix EVD::inverse() const {
    int n = V.sizerow();
    pp::matrix W_inv(n, n);
    for (int i = 0; i < n; i++) {
        if (std::abs(w[i]) < 1e-10) throw std::runtime_error("Matrix is singular!");
        W_inv(i, i) = 1.0 / w[i];
    }
    return V * W_inv * V.transpose(); // A⁻¹ = V W⁻¹ V^T
}
    


}
 //namespace pp
