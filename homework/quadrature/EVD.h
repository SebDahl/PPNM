// #ifndef EVD_H
// #define EVD_H

// #include <iostream>
// #include <vector>
// #include <cmath>
// #include "matrix.h"

// namespace pp {
//     void timesJ(pp::matrix& A, int p, int q, double theta);
//     void Jtimes(pp::matrix& A, int p, int q, double theta);

//     class EVD {
//     private:
//         pp::vector w;
//         pp::matrix V;
//     public:
//         EVD(const pp::matrix& A);  // Constructor (decomposition)
//         pp::vector getW() const;   // Getter for w
//         pp::matrix getV() const;   // Getter for V
//     };
// } // namespace pp

// #endif


#ifndef EVD_H
#define EVD_H

#include <iostream>
#include <vector>
#include <cmath>
#include "matrix.h"

namespace pp {

    class EVD {
    private:
        pp::vector w;  // Eigenvalues
        pp::matrix V;  // Eigenvectors
        pp::matrix D;  // Diagonal eigenvalue matrix

        // Helper functions for Jacobi rotations
        static void timesJ(pp::matrix& A, int p, int q, double theta);
        static void Jtimes(pp::matrix& A, int p, int q, double theta);
        void cyclic_jacobi(pp::matrix& A);

    public:
        // Constructor (performs diagonalization)
        EVD(const pp::matrix& A);

        // Getters
        pp::vector getW() const { return w; }
        pp::matrix getV() const { return V; }
        pp::matrix getD() const { return D; }

        // Additional methods
        pp::matrix inverse() const;
        double det() const;
    };

} // namespace pp

#endif
