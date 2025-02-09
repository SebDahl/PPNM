#include <iostream>
#include <limits>
#include <string>

int main() {
    std::cout << "Maximum value for int" << std::endl;
    int i = 1;
    while (i + 1 > i && i > 0) { // Stop when i + 2 would overflow
        i++;
    }
    std::cout << i-1  << std::endl; // Prints the maximum value for int

    // Compare to expected value
    std::cout << (i - 1 == std::numeric_limits<int>::max() ? "pass" : "fail") << std::endl;

    std::cout << "Minimum value for int" << std::endl;
    int j = -1;
    while (j - 1 < j && j < 0) { // Stop when i + 2 would overflow
        j--;
    }
    std::cout << j  +1<< std::endl; // Prints the maximum value for int

    // Compare to expected value
    std::cout << (j +1 == std::numeric_limits<int>::min() ? "pass" : "fail") << std::endl;
    return 0;
}