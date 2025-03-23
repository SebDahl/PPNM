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
#include <stdexcept>


#include <vector>

struct qspline {
    pp::vector x, y, b, c;

    // Constructor
    qspline(const pp::vector& xs, const pp::vector& ys)
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
int binsearch(pp::vector x, double z)
	{/* locates the interval for z by bisection */ 
	if( z<x[0] || z>x[x.size()-1] ){ throw std::runtime_error("binsearch: bad z");}
	int i=0, j=x.size()-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
	}
double linterp(pp::vector x, pp::vector y, double z){
    int i=binsearch(x,z);
    double dx=x[i+1]-x[i]; if(!(dx>0)){ throw std::runtime_error("uups...");}
    double dy=y[i+1]-y[i];
    return y[i]+dy/dx*(z-x[i]);
    }



double linterpIntegrate(pp::vector x, pp::vector y, double z) {
    if (z < x[0] || z > x[x.size() - 1]) {
        throw std::runtime_error("linterpIntegrate: z out of bounds");
    }
    
    double integral = 0.0;
    for(int i=0;x[i+1] < z; i++;) {
        double dx = x[i+1] - x[i];
        double dy = y[i+1] - y[i];
        double slope = dy / dx;
        integral += y[i] * dx + 0.5 * slope * dx * dx;
    }
    
    double dx = z - x[i];
    double dy = y[i+1] - y[i];
    double slope = dy / (x[i+1] - x[i]);
    integral += y[i] * dx + 0.5 * slope * dx * dx;
    
    return integral;
}




int main(){

    pp::vector x_i(10);
    pp::vector y_i(10);
    for (int i = 0; i < 10; i++) {
        x_i[i] = i;
        y_i[i] = std::cos(x_i[i]);
    }

    pp::vector::write(x_i, "x_i.txt");
    pp::vector::write(y_i, "y_i.txt");

    double z = 3.5;
    double y = linterp(x_i, y_i, z);
    std::cout << "Interpolated value at z = " << z << " is " << y << std::endl;
    pp::vector x_interpolated(0);
    pp::vector x_integrated(0);
    pp::vector y_interpolated(0);
    pp::vector y_integrated(0);
    std::cout << "Test" << x_i[9] << std::endl;

    for(double i = 0; i < 9; i += 0.1){
        x_interpolated.append(i);
        x_integrated.append(i);
        y_interpolated.append(linterp(x_i, y_i, i));
        y_integrated.append(linterpIntegrate(x_i, y_i, i));
    }
    pp::vector::write(x_interpolated, "x_interpolated.txt");
    pp::vector::write(y_interpolated, "y_interpolated.txt");
    pp::vector::write(x_integrated, "x_integrated.txt");
    pp::vector::write(y_integrated, "y_integrated.txt");
    

    return 0;
}








