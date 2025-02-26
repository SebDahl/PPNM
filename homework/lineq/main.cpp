#include <iostream>
#include <vector>
#include <cmath>
#include "matrix.h"
#include "QR.h"
#include <random>

int n = 20;
int m = 10;

int main() {
    pp::matrix A = pp::matrix(m,n);
    std::random_device rd;  
    std::uniform_int_distribution<int> dist(1, 10); 

    for (int i=0; i<n; i++){
        A(i,i) =  



    }




    return 0;
}
