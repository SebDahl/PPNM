#include <iostream>
#include <cmath>
#include "sfuns.h"
int main(int argc, char* argv[]){
    std::cout << "Hello" << std::endl;
    double x =std::sin(1.0);
    double y = std::cos(1.0);
    double z = sfuns::fgamma(1.0);
    std::cout << "x=" << x << " y=" << y << " z=" << z<< std::endl;
    return 0;
}