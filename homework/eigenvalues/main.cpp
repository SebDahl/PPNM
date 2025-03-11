#include <iostream>
#include <cmath>
#include "matrix.h"
#include "EVD.h"
#include "QR.h"
#include <random>
#include <fstream>
#include <string>




int main(){
    int n = 10;
    int m = n/2;

    // Generate random symmetric matrix A
    pp::matrix A(n, n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(1.0, 10.0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double rand = dist(gen);
            A(i, j) = rand;
            A(j, i) = rand;
        }
    }

    // Perform cyclic sweep
    pp::EVD EVD(A);
    pp::vector w = EVD.getW();
    pp::matrix V = EVD.getV();
    pp::matrix D = EVD.getD();

    // write w to "w.txt" using the function write from matrix.h
    pp::vector::write(w, "w.txt");

    // write V to "V.txt" using the function write from matrix.h
    pp::matrix::write(V, "V.txt");

    // write D to "D.txt" using the function write from matrix.h
    
    pp::matrix::write(D, "D.txt");


    // Verification
    pp::matrix VT = V.transpose();




    pp::matrix VTA = VT * A;
    pp::matrix VTAV = VTA * V;

    std::cout << "Checking if V^T * A * V = D: " << (pp::approx_equal(VTAV, D) ? "PASSED" : "FAILED") << std::endl;

    pp::matrix::write(VTAV, "VTAV.txt");

    pp::matrix VD = V * D;
    pp::matrix VDVT = VD * VT;

    std::cout << "Checking if V * D * V^T = A: " << (pp::approx_equal(VDVT, A) ? "PASSED" : "FAILED") << std::endl;


    // check that V^T * V = I
    pp::matrix VTV = VT * V;
    bool isOrthogonal = true;
    for (int i = 0; i < VTV.sizerow(); i++) {
        for (int j = 0; j < VTV.sizecol(); j++) {
            if (i == j && std::abs(VTV(i, j) - 1.0) > 1e-10) {
                isOrthogonal = false;
                break;
            }
            if (i != j && std::abs(VTV(i, j)) > 1e-10) {
                isOrthogonal = false;
                break;
            }
        }
    }

    std::cout << "Checking if V^T * V = I: " << (isOrthogonal ? "PASSED" : "FAILED") << std::endl;

    pp::matrix VVT = V * VT;
    bool isOrthogonal2 = true;
    for (int i = 0; i < VVT.sizerow(); i++) {
        for (int j = 0; j < VVT.sizecol(); j++) {
            if (i == j && std::abs(VVT(i, j) - 1.0) > 1e-10) {
                isOrthogonal2 = false;
                break;
            }
            if (i != j && std::abs(VVT(i, j)) > 1e-10) {
                isOrthogonal2 = false;
                break;
            }
        }
    }
    std::cout << "Checking if V * V^T = I: " << (isOrthogonal2 ? "PASSED" : "FAILED") << std::endl;

    return 0;
}








