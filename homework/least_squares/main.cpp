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





pp::vector lsfit(std::vector<std::function<double(double)>> fs, 
                 pp::vector x, pp::vector y, pp::vector dy) {
    int n = x.size();
    int m = fs.size();

    pp::matrix A(n, m);
    pp::vector b(n);
    pp::matrix W(n, n); // Diagonal weight matrix

    // Fill A and b, and set W as 1/dy^2
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            A(i, j) = fs[j](x[i]) / dy[i]; // Apply weighting
        }
        b[i] = y[i] / dy[i]; // Apply weighting
        W(i, i) = 1.0 / (dy[i] * dy[i]); // Set weight matrix
    }

    // Perform QR decomposition on A
    pp::QR qrA(A);
    pp::vector c = qrA.solve(b); // Solve for c

    return c;
}
std::pair<pp::vector, pp::vector> lsfit_with_uncertainty(std::vector<std::function<double(double)>> fs, 
                                                          pp::vector x, pp::vector y, pp::vector dy) {
    int n = x.size();
    int m = fs.size();

    pp::matrix A(n, m);
    pp::vector b(n);
    pp::matrix W(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            A(i, j) = fs[j](x[i]) / dy[i];  // Weight each function
        }
        b[i] = y[i] / dy[i];  // Weight the observations
        W(i, i) = 1.0 / (dy[i] * dy[i]);  // Diagonal weight matrix
    }

    // Compute A^T W A and A^T W b
    pp::matrix ATW = A.transpose() * W;
    pp::matrix ATA = ATW * A;
    pp::vector ATb = ATW * b;

    // Solve (A^T W A) c = A^T W b using QR decomposition
    pp::QR qrATA(ATA);
    pp::vector c = qrATA.solve(ATb);

    // Solve for covariance matrix: (A^T W A) C = I
    pp::matrix I = pp::matrix::identity(m);
    pp::matrix C(m, m);
    for (int j = 0; j < m; j++) {
        pp::vector col = qrATA.solve(I.get_col(j));  // Solve for each column of C
        C.set_col(j, col);
    }

    // Compute uncertainties in fit parameters: sigma_c[j] = sqrt(C[j,j])
    pp::vector sigma_c(m);
    for (int j = 0; j < m; j++) {
        sigma_c[j] = std::sqrt(C(j, j));
    }

    return {c, sigma_c}; // Return coefficients and their uncertainties
}

// std::vector lsfit(std::function<double, double>[] fs, vector x, vector y, vector dy, pp::vector c_2);
double plot_func(double z, pp::vector c_2) {
    return std::exp(c_2[0]) * std::exp(-c_2[1]*z);
};

int main(){

    std::vector<std::function<double(double)>> fs = {
        [](double z) { return 1.0; },
        [](double z) { return z; },
        [](double z) { return z * z; }
        };
    
    std::vector<std::function<double(double)>> fs_exp = {
        [](double z) { return 1.0; },
        [](double z) { return -z; },
        // [](double z) { return z * z; }
        };

    pp::matrix Data = pp::matrix::loadtxt("data.data");
    pp::vector x = Data.get_col(0);
    pp::vector y = Data.get_col(1);
    pp::vector dy = Data.get_col(2);
    

    pp::matrix A = pp::matrix(x.size(), fs.size());
    pp::vector b = y;
    // std::cout << "Size of x: " << x.size() << std::endl;

    pp::matrix::write(A, "A.txt");
    pp::vector::write(x, "x.txt");
    // std::cout << "x.size() = " << x.size() << ", fs.size() = " << fs.size() << std::endl;

    
    for (int i = 0; i < fs.size(); i++) {
    for (int j = 0; j < x.size(); j++) {
        double value = fs[i](x[j]);
        // std::cout << "Setting A(" << j << ", " << i << ") to " << value << std::endl;
        A(j, i) = value;
    }
    }
    pp::matrix::write(A, "A2.txt");

    pp::QR qrA(A);
    pp::vector c = qrA.solve(b);

    pp::vector yfit = A.transpose() * c;
    pp::vector residuals = y - yfit;
    pp::vector chi2 = residuals * residuals / (dy * dy);

    pp::vector::write(c, "c.txt");

    // std::cout << "Coefficients: " << c << std::endl;
    // std::cout << "Residuals: " << residuals << std::endl;
    // std::cout << "Chi^2: " << chi2[0] << std::endl;

    pp::matrix data_rad = pp::matrix::loadtxt("data_rad.data");
    pp::vector x_rad = data_rad.get_col(0);

    pp::vector y_rad = data_rad.get_col(1);
    pp::vector dy_rad = data_rad.get_col(2);


    for (int i=0; i < x_rad.size(); i++) {
        y_rad[i] = std::log(y_rad[i]);
        dy_rad[i] = dy_rad[i] / y_rad[i];
    }

    pp::vector c_2 = lsfit(fs_exp, x_rad, y_rad, dy_rad);
    pp::vector::write(c_2, "c_2.txt");

    pp::vector plot_data = pp::vector(x_rad.size());
    for (int i = 0; i < x_rad.size(); i++) {
        plot_data[i] = plot_func(x_rad[i], c_2);
    }
    pp::vector::write(plot_data, "plot_data.txt");

    pp::vector c_3 = lsfit_with_uncertainty(fs_exp, x_rad, y_rad, dy_rad).first;
    pp::vector dc_3 = lsfit_with_uncertainty(fs_exp, x_rad, y_rad, dy_rad).second;
    pp::vector::write(c_3, "c_3.txt");
    pp::vector::write(dc_3, "dc_3.txt");
    std::cout << "The decay constant is " << c_3[1] << " +/- " << dc_3[1] << " days^-1." << std::endl;


    std::cout << "The half-life of the isotope is " << std::log(2)/c_3[1] << " +/- " << std::log(2)/c_3[1] * dc_3[1] << " days." << std::endl;
    std::cout << "The modern value of the half-life is 3.6316(14) days, so it does not agree within uncertainty but not too far." << std::endl;


    pp::vector c_3_upper = c_3 + dc_3;
    pp::vector c_3_lower = c_3 - dc_3;

    pp::vector plot_data_upper = pp::vector(x_rad.size());
    pp::vector plot_data_lower = pp::vector(x_rad.size());

    for (int i = 0; i < x_rad.size(); i++) {
        plot_data_upper[i] = plot_func(x_rad[i], c_3_upper);
        plot_data_lower[i] = plot_func(x_rad[i], c_3_lower);
    }

    pp::vector::write(plot_data_upper, "plot_data_upper.txt");
    pp::vector::write(plot_data_lower, "plot_data_lower.txt");

    

    

    return 0;
}








