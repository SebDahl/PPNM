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


static pp::vector plainmc(std::function<double(pp::vector)> f, pp::vector a, pp::vector b, int N){
    int dim = a.size();
    double V=1;
    for(int i=0;i<dim;i++){V*=b[i]-a[i];}
    double sum=0, sum2=0;
    auto x = pp::vector(dim);
    auto rnd = std::

}


int main(){



    return 0;
}








