#include <iostream>
#include <cmath>
#include "matrix.h"
#include "sfuns.h"
#include "EVD.h"
#include "QR.h"
#include <random>
#include <fstream>
#include <functional>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>




static pp::matrix jacobian(
    std::function<pp::vector(pp::vector)> f,
    pp::vector x,
    pp::vector fx,
    pp::vector deltax = 0)
    {
        if (deltax.size() == 0) {
            deltax.resize(x.size());
            for (int i = 0; i < x.size(); ++i) {
                deltax[i] = std::abs(x[i]) * std::pow(2.0, -26);
            }
        }
    if(fx.size() != x.size()){
        fx = f(x);
    }
    pp::matrix J(x.size(), x.size());
    for (int j=0; j < x.size(); j++){
        x[j] += deltax[j];
        pp::vector df = f(x)- fx;
        for (int i=0; i < x.size(); i++){
            J[i,j] = df[i]/deltax[j];
        }
        x[j] -= deltax[j];
    }
    return J;
    }

static pp::vector newtons_method_back(
    std::function<pp::vector(pp::vector)> f,
    pp::vector start, 
    double accuracy = 1e-2,
    double lambda_min = 1e-3,
    pp::vector deltax = pp::vector{0}){
        pp::vector x = start;
        pp::vector fx = f(x), z, fz;
        do{/*Newtons iterations*/
            if(fx.norm() < accuracy) break;
            pp::matrix J = jacobian();
            auto QRJ = pp::QR(J);
            pp::vector Dx = QRJ.solve(-1*fx);
            double lambda = 1;
            do{/*line search*/
                z = x + lambda*Dx;
                fz = f(z);
                if(fz.norm() < (1 - lambda/2)*fx.norm() ) break;
                if(lambda < lambda_min) break;
                lambda *= 0.5;
            }while(true);
        x=z; fx=fz;
        }while(true);
    return x;

}


int main(){



    return 0;
}



