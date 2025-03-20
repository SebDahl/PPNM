#ifndef QR_H
#define QR_H

#include <iostream>
#include <vector>
#include <cmath>
#include "matrix.h"
namespace pp {

class QR {
private:
    pp::matrix Q;
    pp::matrix R;


public:
    QR(const pp::matrix& A); // Constructor (decomp)
    pp::matrix getQ() const; // Getter for Q
    pp::matrix getR() const; // Getter for R
    pp::vector solve(const pp::vector& b) const; 
    pp::matrix inverse(const pp::matrix& A) const; 
    double det() const;
    

};

} //namespace pp

#endif


