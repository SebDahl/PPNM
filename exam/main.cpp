// main.cpp (modified to include least-squares missing samples recovery)
#include <iostream>
#include <cmath>
#include "matrix.h"
#include "EVD.h"
#include "QR.h"
#include <random>
#include <fstream>
#include <functional>
#include <string>
#include <cmath>

// Construct finite-difference matrix D for second derivatives
// pp::matrix D_matrix(int n) {
//     pp::matrix D(n - 2, n);
//     for (int i = 0; i < n - 2; i++) {
//         D(i, i) = 1.0;
//         D(i, i + 1) = -2.0;
//         D(i, i + 2) = 1.0;
//     }
//     return D;
// }
pp::matrix D_matrix(const pp::vector& x){
    int n = x.size();
    std::cout << "Size of x: " << n << std::endl;
    pp::matrix D(n, n);
    D(0, 0) = 0.0; // First element is 1
    D(0, 1) = -0.0; // Second element is -2
    D(0, 2) = 0.0; // Third element is 1
    D(n-1, n-1) = 0.0; // Last element is 1
    D(n-1, n-2) = -0.0; // Second to last element is -2
    D(n-1, n-3) = 0.0; // Third to last element is 1

    // D(0, 0) = 0.0; // First element is 1
    // D(0, 1) = -0.0; // Second element is -2
    // D(0, 2) = 0.0; // Third element is 1
    // D(n-1, n-1) = 0.0; // Last element is 1
    // D(n-1, n-2) = -0.0; // Second to last element is -2
    // D(n-1, n-3) = 0.0; // Third to last element is 1
    for (int i = 1; i < n-1; i++){
        D(i, i-1) = 1.0;
        D(i, i) = -2.0;
        D(i, i+1) = 1.0;
    }
    return D;
}

pp::matrix D3_matrix(const pp::vector& x){
    int n = x.size();
    pp::matrix D(n, n);
    for (int i = 0; i < 2; i++) {
    D(i, 0) = -1.0; // First element is 1
    D(i, 1) = 3.0; // Second element is -2
    D(i, 2) = -3.0; // Third element is 3
    D(i, 2) = 1.0; // Third element is 1    
    }
    for (int i = n-3; i < n; i++) {
    D(i, n-1) = 1.0; // Last element is 1
    D(i, n-2) = -3.0; // Second to last element is -2
    D(i, n-3) = 3.0; // Third to last element is 3
    D(i, n-4) = -1.0; // Fourth to last element is 1
    }

    for (int i = 2; i < n-2; i++){
        D(i, i-2) = -0.5;
        D(i, i-1) = 1;
        D(i, i) = 0.0;
        D(i, i+1) = -1.0;
        D(i, i+2) = 0.5;
    }
    return D;
}



// Construct insertion matrix M for missing indices
pp::matrix M_matrix(int N, const std::vector<int>& missing_indices) {
    int m = missing_indices.size();
    pp::matrix M(N, m);
    for (int k = 0; k < m; ++k) {
        int idx = missing_indices[k];
        M(idx, k) = 1.0;
    }
    return M;
}

// Perform least-squares missing sample recovery
std::pair<pp::vector, pp::vector> recover_missing(pp::vector y) {
    int N = y.size();
    std::vector<int> missing_indices = {};
    for (int i = 0; i < N; i++){
        if (y[i] == 0.0) { // missing samples are marked as zero
            missing_indices.push_back(i);
        }
    }
    int m = missing_indices.size();

    // Build M and D
    pp::matrix M = M_matrix(N, missing_indices);
    pp::matrix D = D_matrix(N);

    // Compute DM and Dy
    pp::matrix DM = D * M;
    pp::vector Dy = D * y;

    // Solve DM z = -Dy
    pp::vector rhs = Dy * (-1.0);
    pp::QR qr(DM);
    pp::vector z = qr.solve(rhs);

    // Recover full signal x = y + M z
    pp::vector x = y;
    for (int k = 0; k < m; ++k) {
        int idx = missing_indices[k];
        x[idx] += z[k];
    }

        // Estimate covariance matrix of z
    pp::matrix DTDM = DM.transpose() * DM;
    pp::QR qr_cov(DTDM);
    pp::matrix C_inv(m, m);
    pp::matrix I = pp::matrix::identity(m);
    for (int j = 0; j < m; ++j) {
        pp::vector col = qr_cov.solve(I.get_col(j));
        C_inv.set_col(j, col);
    }

    // Estimate residual variance sigma^2
    pp::vector residual = DM * z + Dy;
    double sigma2 = (residual * residual) / (DM.sizerow() - m);

    // Final covariance matrix and uncertainties
    pp::matrix Cov_z = C_inv * sigma2;
    pp::vector sigma_z(m);
    for (int i = 0; i < m; ++i)
        sigma_z[i] = std::sqrt(Cov_z(i, i));

    return {x, sigma_z};
}


// pp::vector declipping(pp::vector y) {
//     int N = y.size();
//     std::vector<int> missing_indices = {};
//     for (int i = 0; i < N; i++){
//         if (y[i] == 0.0) { // missing samples are marked as zero
//             missing_indices.push_back(i);
//         }
//     }
//     int m = missing_indices.size();

//     // Build M and D
//     pp::matrix M = M_matrix(N, missing_indices);
//     pp::matrix D = D3_matrix(N);

//     // Compute DM and Dy
//     pp::matrix DM = D * M;
//     pp::vector Dy = D * y;

//     // Solve DM z = -Dy
//     pp::vector rhs = Dy * (-1.0);
//     pp::QR qr(DM);
//     pp::vector z = qr.solve(rhs);

//     // Recover full signal x = y + M z
//     pp::vector x = y;
//     for (int k = 0; k < m; ++k) {
//         int idx = missing_indices[k];
//         x[idx] += z[k];
//     }

//     return x;
// }

int main() {
    // Create test signal with missing samples
    int N = 20;
    pp::vector y(N);
    for (int i = 0; i < N; ++i)
        y[i] = std::sin(0.3 * i)*i;
    pp::vector::write(y, "actual_signal.txt");

    std::vector<int> missing = {0,1,2, 5, 6, 7, 8, 19};
    for (int i : missing)
        y[i] = 0.0;  // Mark missing samples as zero

    pp::vector x_recovered = recover_missing(y).first;
    pp::vector sigma_z = recover_missing(y).second;
    pp::vector::write(sigma_z, "uncertainties.txt");
    pp::vector missing_indices(missing.size());
    for (int i = 0; i < missing.size(); ++i) {
        missing_indices[i] = missing[i];
    }
    pp::vector::write(missing_indices, "missing_indices.txt");


    pp::vector::write(y, "signal_with_missing.txt");
    pp::vector::write(x_recovered, "recovered_signal.txt");

    std::cout << "Recovered signal written to 'recovered_signal.txt'\n";

    std::ofstream errfile("errorbars.txt");
    for (int i = 0; i < missing.size(); ++i) {
        int idx = missing[i];
        errfile << idx << " " << x_recovered[idx] << " " << sigma_z[i] << "\n";
    }
    errfile.close();

    // pp::vector y2(N);
    // for (int i = 0; i < N; ++i)
    //     y2[i] = std::sin(2 * i);
    // pp::vector declipped_signal = declipping(y2);
    // pp::vector::write(y2, "signal_with_clipping.txt");
    // pp::vector::write(declipped_signal, "declipped_signal.txt");

    
    return 0;
}
