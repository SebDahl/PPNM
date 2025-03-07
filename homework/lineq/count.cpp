#include <iostream>
#include <random>
#include <chrono>
#include "matrix.h"
#include "QR.h"
#include <cstdlib> // For atoi

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <matrix_size>\n";
        return 1;
    }
    int n = std::atoi(argv[1]);
    if(n <= 0) {
        std::cerr << "Matrix size must be positive\n";
        return 1;
    }

    // Generate random matrix A of size n x n
    pp::matrix A(n, n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(1.0, 10.0);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            A(i, j) = dist(gen);
        }
    }

    // Time the QR decomposition (if you want internal timing, you could also use chrono here)
    pp::QR qr(A);

    return 0;
}
