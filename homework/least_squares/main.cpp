#include <iostream>
#include <cmath>
#include "matrix.h"
#include "EVD.h"
#include "QR.h"
#include <random>
#include <fstream>
#include <functional>
#include <string>





pp::vector lsfit(std::vector<std::function<double(double)>> fs_2, vector x, vector y, vector dy) {
    pp::matrix A = pp::matrix(x.size(), fs_2.size());
    pp::vector b = y;

    for (int i = 0; i < fs.size(); i++) {
        for (int j = 0; j < x.size(); j++) {
            double value = fs[i](x[j]);
            // std::cout << "Setting A(" << j << ", " << i << ") to " << value << std::endl;
            A(j, i) = value;
        }
        }

    pp::QR qrA(A);
    pp::vector c = qrA.solve(b);

    return c;
}

// std::vector lsfit(std::function<double, double>[] fs, vector x, vector y, vector dy);


int main(){

    std::vector<std::function<double(double)>> fs = {
        [](double z) { return 1.0; },
        [](double z) { return z; },
        [](double z) { return z * z; }
        };
    
    std::vector<std::function<double(double)>> fs_exp = {
        [](double z) { return 1.0; },
        [](double z) { return z; },
        [](double z) { return z * z; }
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

    
    // pp::vector c_2 = lsfit(fs, x, y, dy);
    // pp::vector::write(c_2, "c_2.txt");

    

    

    return 0;
}








