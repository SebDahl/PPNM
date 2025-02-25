#ifndef QR_H
#define QR_H

#include <iostream>
#include <vector>
#include <cmath>
#include "matrix.h"

class QR {
private:
    std::matrix<double> Q;
    std::matrix<double> R;


public:
    QR(const std::matrix<double>& A); // Constructor (decomp)
    std::matrix<double> getQ() const; // Getter for Q
    std::matrix<double> getR() const; // Getter for R
    std::matrix<double> solve(const std::matrix<double>& b) const; 
    double det() const;
    

};

#endif


