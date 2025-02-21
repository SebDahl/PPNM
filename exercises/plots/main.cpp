#include <iostream>
#include <cmath>
#include "sfuns.h"

int main() {
    for (double x = -2; x <= 2; x += 0.1) {
        std::cout << x << " " << sfuns::erf(x) << std::endl;
    }
    return 0;
}
