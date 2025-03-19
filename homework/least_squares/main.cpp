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





pp::vector lsfit(std::vector<std::function<double(double)>> fs_2, pp::vector x_2, pp::vector y_2, pp::vector dy_2) {
    pp::matrix A_2 = pp::matrix(x_2.size(), fs_2.size());
    pp::vector b_2 = y_2;

    for (int i = 0; i < fs_2.size(); i++) {
        for (int j = 0; j < x_2.size(); j++) {
            double value = fs_2[i](x_2[j]);
            // std::cout << "Setting A(" << j << ", " << i << ") to " << value << std::endl;
            A_2(j, i) = value;
        }
        }

    pp::QR qrA(A_2);
    pp::vector c = qrA.solve(b_2);

    return c;
}

// std::vector lsfit(std::function<double, double>[] fs, vector x, vector y, vector dy);
// double plot_func(double z){
//     return std::exp(c_2[0] + c_2[1] * z);
// };

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
    std::cout << "Size of x: " << x.size() << std::endl;

    pp::matrix::write(A, "A.txt");
    pp::vector::write(x, "x.txt");
    std::cout << "x.size() = " << x.size() << ", fs.size() = " << fs.size() << std::endl;

    
    for (int i = 0; i < fs.size(); i++) {
    for (int j = 0; j < x.size(); j++) {
        double value = fs[i](x[j]);
        // std::cout << "Setting A(" << j << ", " << i << ") to " << value << std::endl;
        A(j, i) = value;
    }
    }
    std::cout << "Test" << std::endl;
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
        x_rad[i] = std::log(x_rad[i]);
        y_rad[i] = std::log(y_rad[i]);
        dy_rad[i] = dy_rad[i] / y_rad[i];
    }

    pp::vector c_2 = lsfit(fs_exp, x_rad, y_rad, dy_rad);
    pp::vector::write(c_2, "c_2.txt");

    // pp::vector plot_data = pp::vector(x_rad.size());


    

    

    return 0;
}








