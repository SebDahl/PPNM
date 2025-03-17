#include <iostream>
#include <cstdlib>  // For atof
#include <string>
#include "matrix.h" // 
#include "EVD.h"    //
#include "QR.h"     //

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
    // std::cout << "Hamiltonian Matrix (H):\n";
    // H.print("Hamiltonian Matrix:", stdout);


    // Now diagonalize the Hamiltonian matrix using EVD
    pp::EVD evd(H);
    pp::vector w = evd.getW();
    pp::matrix V = evd.getV();
    pp::matrix D = evd.getD();

    pp::matrix::write(D, "H_D.txt");
    pp::vector::write(w, "H_w.txt");
    pp::matrix::write(V, "H_V.txt");
    std::cout << "finished writing files, now moving on to convergence" << std::endl;

    double E0 = w[0];
    pp::vector E0s(0);
    pp::vector rmaxs(0);
    for (int j = 5; j < 15; j++) {
        rmax = j;
        dr = 0.3;

    
        npoints = static_cast<int>(rmax / dr) - 1;
        pp::vector r(npoints);
    
        for (int i = 0; i < npoints; i++) {
            r[i] = dr * (i + 1);
        }
    
        // Create Hamiltonian matrix using pp::matrix
        pp::matrix H(npoints, npoints);
        
        factor = -0.5 / (dr * dr);
        
        for (int i = 0; i < npoints - 1; i++) {
            H(i, i)     = -2 * factor;
            H(i, i + 1) = factor;
            H(i + 1, i) = factor;
        }
        H(npoints - 1, npoints - 1) = -2 * factor;
    
        for (int i = 0; i < npoints; i++) {
            H(i, i) += -1.0 / r[i];
        }
    
    
    
        // Now diagonalize the Hamiltonian matrix using EVD
        pp::EVD evd(H);
        pp::vector w = evd.getW();
        // pp::matrix V = evd.getV();
        // pp::matrix D = evd.getD();
    
        E0 = w[0];
        
        E0s.append(E0);
        rmaxs.append(rmax);


    }
    pp::vector::write(E0s, "E0s_rmax.txt");
    pp::vector::write(rmaxs, "rmaxs.txt");

    pp::vector E0s2(0);
    pp::vector drs(0);
    for (double k = 0.1; k < 0.7; k += 0.1) {
        rmax = 10;
        dr = k;

    
        npoints = static_cast<int>(rmax / dr) - 1;
        pp::vector r(npoints);
    
        for (int i = 0; i < npoints; i++) {
            r[i] = dr * (i + 1);
        }
    
        // Create Hamiltonian matrix using pp::matrix
        pp::matrix H(npoints, npoints);
        
         factor = -0.5 / (dr * dr);
        
        for (int i = 0; i < npoints - 1; i++) {
            H(i, i)     = -2 * factor;
            H(i, i + 1) = factor;
            H(i + 1, i) = factor;
        }
        H(npoints - 1, npoints - 1) = -2 * factor;
    
        for (int i = 0; i < npoints; i++) {
            H(i, i) += -1.0 / r[i];
        }
    
    
    
        // Now diagonalize the Hamiltonian matrix using EVD
        pp::EVD evd(H);
        pp::vector w = evd.getW();
        // pp::matrix V = evd.getV();
        // pp::matrix D = evd.getD();
    
        E0 = w[0];

        E0s2.append(E0);
        drs.append(dr);

    }
    pp::vector::write(E0s2, "E0s_dr.txt");
    pp::vector::write(drs, "drs.txt");

    std::cout << "Convergence test complete" << std::endl;



    return 0;
}
