#include <iostream>
#include <cstdlib>  // For atof
#include <string>
#include "matrix.h" // 

int main(int argc, char* argv[]) {
    double rmax = 0.0, dr = 0.0;

    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-rmax" && i + 1 < argc) {
            rmax = std::atof(argv[i + 1]);
            i++;
        } else if (arg == "-dr" && i + 1 < argc) {
            dr = std::atof(argv[i + 1]);
            i++;
        }
    }

    // Validate input
    if (rmax <= 0 || dr <= 0) {
        std::cerr << "Usage: " << argv[0] << " -rmax <value> -dr <value>\n";
        return 1;
    }

    int npoints = static_cast<int>(rmax / dr) - 1;
    pp::vector r(npoints);

    for (int i = 0; i < npoints; i++) {
        r[i] = dr * (i + 1);
    }

    // Create Hamiltonian matrix using pp::matrix
    pp::matrix H(npoints, npoints);
    
    double factor = -0.5 / (dr * dr);
    
    for (int i = 0; i < npoints - 1; i++) {
        H(i, i)     = -2 * factor;
        H(i, i + 1) = factor;
        H(i + 1, i) = factor;
    }
    H(npoints - 1, npoints - 1) = -2 * factor;

    for (int i = 0; i < npoints; i++) {
        H(i, i) += -1.0 / r[i];
    }

    // Print the Hamiltonian matrix for verification
    std::cout << "Hamiltonian Matrix (H):\n";
    H.print("Hamiltonian Matrix:", stdout);

    return 0;
}
