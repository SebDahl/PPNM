#include <iostream>
#include <cmath>
#include "matrix.h"
#include "QR.h"
#include <random>
#include <fstream>
#include <string>

bool isUpperTriangular(const pp::matrix& R) {
    int rows = R.sizerow();
    int cols = R.sizecol();
    for (int i = 1; i < rows; i++) {  // Start from row 1
        for (int j = 0; j < i; j++) { // Check below diagonal
            if (std::abs(R(i, j)) > 1e-10) {
                return false;  // Non-zero entry below diagonal
            }
        }
    }
    return true;
}

bool isOrthogonal(const pp::matrix& Q) {
    int n = Q.sizerow();
    int m = Q.sizecol();
    pp::matrix QT = Q.transpose();
    pp::matrix QTQ = QT * Q;

    // Check if QTQ is approximately the identity matrix
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (i == j && std::abs(QTQ(i, j) - 1.0) > 1e-10) return false; // Check diagonal
            if (i != j && std::abs(QTQ(i, j)) > 1e-10) return false; // Check off-diagonal
        }
    }
    return true;
}

bool isQRReconstructionCorrect(const pp::matrix& Q, const pp::matrix& R, const pp::matrix& A) {
    pp::matrix QR = Q * R;
    int rows = A.sizerow();
    int cols = A.sizecol();

    // Check if QR is approximately equal to A
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (std::abs(QR(i, j) - A(i, j)) > 1e-10) {
                return false;
            }
        }
    }
    return true;
}

int main() {
    int m = 20;  // Rows (tall matrix)
    int n = 10;  // Columns (m > n)

    // Generate random matrix A
    pp::matrix A(m, n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(1.0, 10.0);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            A(i, j) = dist(gen);
        }
    }

    std::cout << "First we perform and check the QR decomposition:" << std::endl;
    // Perform QR decomposition
    pp::QR qr(A);
    pp::matrix Q = qr.getQ();
    pp::matrix R = qr.getR();

    // Validate the decomposition
    std::cout << "Checking if R is upper triangular: " << (isUpperTriangular(R) ? "PASSED" : "FAILED") << std::endl;
    std::cout << "Checking if Q^T * Q = I: " << (isOrthogonal(Q) ? "PASSED" : "FAILED") << std::endl;
    std::cout << "Checking if QR = A: " << (isQRReconstructionCorrect(Q, R, A) ? "PASSED" : "FAILED") << std::endl;
    // Write Q to "Q.txt"
    std::ofstream qfile("Q.txt");
    if (qfile.is_open()) {
        for (int i = 0; i < Q.sizerow(); i++) {
            for (int j = 0; j < Q.sizecol(); j++) {
                qfile << Q(i, j) << " ";
            }
            qfile << "\n";
        }
        qfile.close();
    } else {
        std::cerr << "Error opening Q.txt\n";
    }


    // Write R to "R.txt"
    std::ofstream rfile("R.txt");
    if (rfile.is_open()) {
        for (int i = 0; i < R.sizerow(); i++) {
            for (int j = 0; j < R.sizecol(); j++) {
                rfile << R(i, j) << " ";
            }
            rfile << "\n";
        }
        rfile.close();
    } else {
        std::cerr << "Error opening R.txt\n";
    }

    // Write A to "A.txt"
    std::ofstream afile("A.txt");
    if (afile.is_open()) {
        for (int i = 0; i < A.sizerow(); i++) {
            for (int j = 0; j < A.sizecol(); j++) {
                afile << A(i, j) << " ";
            }
            afile << "\n";
        }
        afile.close();
    } else {
        std::cerr << "Error opening A.txt\n";
    }



    std::cout << "QR decomposition completed. Now moving on to solving a system of equations, using a square matrix B." << std::endl; 
    // Generate random matrix B
    pp::matrix B(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B(i, j) = dist(gen);
        }
    }

    // Random vector b of same size: n
    pp::vector b(n);
    for (int i = 0; i < n; i++) {
        b[i] = dist(gen);
    }

    // QR decomp of B
    pp::QR qrB(B);
    pp::matrix QB = qrB.getQ();
    pp::matrix RB = qrB.getR();

    // Solve system of equations
    pp::vector x = qrB.solve(b);

    // Validate the solution
    pp::vector Bx = B * x;
    std::cout << "Checking if Bx = b: " << (Bx == b ? "PASSED" : "FAILED") << std::endl;

    // Write B to "B.txt"
    std::ofstream bfile("B.txt");
    if (bfile.is_open()) {
        for (int i = 0; i < B.sizerow(); i++) {
            for (int j = 0; j < B.sizecol(); j++) {
                bfile << B(i, j) << " ";
            }
            bfile << "\n";
        }
        bfile.close();
    } else {
        std::cerr << "Error opening B.txt\n";
    }

    // Write b to "b.txt"
    std::ofstream bfile2("b.txt");
    if (bfile2.is_open()) {
        for (int i = 0; i < b.size(); i++) {
            bfile2 << b[i] << " ";
        }
        bfile2.close();
    } else {
        std::cerr << "Error opening b.txt\n";
    }

    // Write x to "x.txt"
    std::ofstream xfile("x.txt");
    if (xfile.is_open()) {
        for (int i = 0; i < x.size(); i++) {
            xfile << x[i] << " ";
        }
        xfile.close();
    } else {
        std::cerr << "Error opening x.txt\n";
    }

    // Write Bx to "Bx.txt"
    std::ofstream bxfile("Bx.txt");
    if (bxfile.is_open()) {
        for (int i = 0; i < Bx.size(); i++) {
            bxfile << Bx[i] << " ";
        }
        bxfile.close();
    } else {
        std::cerr << "Error opening Bx.txt\n";
    }


    return 0;
}
