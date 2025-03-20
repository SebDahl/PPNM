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


#include <vector>

struct qspline {
    std::vector<double> x, y, b, c; // Use std::vector<double>

    // Constructor
    qspline(const std::vector<double>& xs, const std::vector<double>& ys)
        : x(xs), y(ys), b(xs.size()), c(xs.size()) 
    {
        // Code to compute b and c
    }

    // Function to evaluate the spline
    double evaluate(double z) const {
        // Implement spline evaluation using x, y, b, c
        return 0.0; // Placeholder
    }
};




int main(){

    double asdasd =1;
    pp::vector asdasd(3);
    

    return 0;
}








