#ifndef QR_H
#define QR_H

#include <iostream>
#include <vector>
#include <cmath>
#include "matrix.h"

class QR {
private:
    pp::matrix<double> Q;
    pp::matrix<double> R;


public:
    QR(const pp::matrix<double>& A); // Constructor (decomp)
    pp::matrix<double> getQ() const; // Getter for Q
    pp::matrix<double> getR() const; // Getter for R
    pp::matrix<double> solve(const pp::matrix<double>& b) const; 
    double det() const;
    

};

#endif


