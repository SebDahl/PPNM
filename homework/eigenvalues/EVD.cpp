#include "QR.h"
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cassert>

namespace pp {
EVD::EVD(const pp::matrix& A){  //w, V constructor
    int m = A.sizerow();
    int n = A.sizecol();
    assert(m == n && "Matrix must be square (m == n)");

    w = pp::vector(n);  // w is n × 1
    V = pp::matrix(n, n);  // V is n × n (upper triangular)
    pp::matrix A_copy = A;



void timesJ(pp::matrix& A, int p, int q, double theta){
    double c = std::cos(theta), double s = std::sin(theta);
    for (int i = 0; i < A.sizerow(); i++){
        double temp = c * A(i, p) - s * A(i, q);
        A(i, q) = s * A(i, p) + c * A(i, q);
        A(i, p) = temp;
    }
}

}





} //namespace pp