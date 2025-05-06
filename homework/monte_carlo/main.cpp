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
    for(int i=0;i<dim;i++){V*=(b[i]-a[i]);}
    double sum=0, sum2=0;
    auto x = pp::vector(dim);
    std::random_device rd;
    std::mt19937 gen(rd());
    for (int i=0; i<N; i++){
        for (int j=0; j<dim; j++){
            x[j] = a[j] + (b[j]-a[j])*std::generate_canonical<double, 10>(gen);
        }
    double fx = f(x);
    sum += fx;
    sum2 += fx*fx;
    }
    double mean = sum/N, sigma = std::sqrt((sum2/N - mean*mean));
    pp::vector result = {mean*V, sigma*V/sqrt(N)};
    return result;
}
double corput(int N, int base){
    double q = 0, bk = 1.0/base;
    while (N > 0) {
        q += (N % base) * bk;
        N /= base;
        bk /= base;
    }
    return q;
}

std::vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
    67, 71, 73, 79, 83, 89, 97};



static pp::vector quasimc(std::function<double(pp::vector)> f, pp::vector a, pp::vector b, int N){
    int dim = a.size();
    double V=1;
    for(int i=0;i<dim;i++){V*=(b[i]-a[i]);}
    double sum=0, sum2=0;
    auto x = pp::vector(dim);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            x[j] = a[j] + (b[j] - a[j]) * (corput(i+1, primes[j]));
        }
        double fx = f(x);
        sum += fx;
        sum2 += fx*fx;
    }
    double mean = sum/N, sigma = std::sqrt((sum2/N - mean*mean));
    pp::vector result = {mean*V, sigma*V/sqrt(N)};
    return result;
}

pp::matrix sampler(int N, pp::vector a, pp::vector b) {
    int dim = a.size();
    pp::matrix x(N, dim);  // N rows (samples), dim columns (dimensions)
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            double u = std::generate_canonical<double, 10>(gen);
            x(i, j) = a[j] + (b[j] - a[j]) * u;
        }
    }

    return x;
}

// static pp::vector stratifiedmc(std::function<double(pp::vector)> f, pp::vector a, pp::vector b, int N, int nmin) {
//     int dim = a.size();
//     double V = 1.0;
//     for (int i = 0; i < dim; i++) V *= (b[i] - a[i]);

//     if (N < nmin) {
//         std::cout << "N is too small, using plainmc instead" << std::endl;
//         return plainmc(f, a, b, N);
//     }
//     int N_remaining = N - nmin;
//     if (N_remaining < 2 * nmin) {
//         std::cout << "Not enough points for stratification, fallback to plain MC" << std::endl;
//         return plainmc(f, a, b, N);
//     }

//     // Sample nmin points and evaluate f
//     pp::matrix X = sampler(nmin, a, b); // shape: nmin × dim
//     pp::vector fxs(nmin);
//     for (int i = 0; i < nmin; ++i) fxs[i] = f(X[i]);

//     // Find dimension with largest sub-variance
//     int best_dim = 0;
//     double max_subvar = -1.0;
//     for (int d = 0; d < dim; ++d) {
//         double mid = 0.5 * (a[d] + b[d]);

//         // Partition fxs into left/right by dimension
//         pp::vector left_vals, right_vals;
//         for (int i = 0; i < nmin; ++i) {
//             if (X(i, d) <= mid)
//                 left_vals.append(fxs[i]);
//             else
//                 right_vals.append(fxs[i]);
//         }

//         if (left_vals.size() == 0 || right_vals.size() == 0) continue;

//         double muL = left_vals.mean();
//         double muR = right_vals.mean();
//         double subvar = 0.5 * (muL - muR) * (muL - muR);

//         if (subvar > max_subvar) {
//             max_subvar = subvar;
//             best_dim = d;
//         }
//     }

//     // Split domain along best_dim
//     double mid = 0.5 * (a[best_dim] + b[best_dim]);
//     pp::vector left_b = b.copy();
//     left_b[best_dim] = mid;
//     pp::vector right_a = a.copy();
//     right_a[best_dim] = mid;

//     // Reuse left/right fxs to estimate subvariances
//     pp::vector left_vals, right_vals;
//     for (int i = 0; i < nmin; ++i) {
//         if (X(i, best_dim) <= mid)
//             left_vals.append(fxs[i]);
//         else
//             right_vals.append(fxs[i]);
//     }

//     double vL = left_vals.var(); // assume var() uses Bessel's correction
//     double vR = right_vals.var();
//     double total_v = vL + vR;

    

    
//     // int N_left = std::max(nmin, static_cast<int>(N_remaining * (vL / (vL + vR))));
//     // int N_right = std::max(nmin, N_remaining - N_left);
//     int N_left = static_cast<int>(N_remaining * (vL / (vL + vR)));
//     int N_right = N_remaining - N_left;

//     // Ensure at least nmin per side or fallback to plainmc if impossible
//     if (N_left < nmin || N_right < nmin) {
//         std::cout << "Not enough samples for recursion, falling back to plainmc" << std::endl;
//         return plainmc(f, a, b, N);
//     }

//     // Recursive calls

    
//     pp::vector left_result = stratifiedmc(f, a, left_b, N_left, nmin);
//     pp::vector right_result = stratifiedmc(f, right_a, b, N_right, nmin);

//     // Combine results
//     double I_total = left_result[0] + right_result[0];
//     double err_total = std::sqrt(left_result[1] * left_result[1] + right_result[1] * right_result[1]);

//     return {I_total, err_total};
// }

static pp::vector stratifiedmc(std::function<double(pp::vector)> f, pp::vector a, pp::vector b, int N, int nmin, int max_depth = 10) {
    // Debug output to track recursion
    static int recursion_count = 0;
    recursion_count++;
    
    if (recursion_count % 100 == 0) {
        std::cout << "Recursion depth: " << recursion_count << ", N=" << N << std::endl;
    }
    int dim = a.size();
    double V = 1.0;
    for (int i = 0; i < dim; i++) V *= (b[i] - a[i]);

    // Early exit conditions
    if (N < nmin || max_depth <= 0) {
        // Not enough samples or reached max depth, use plainmc
        std::cout << "Early exit: N=" << N << ", nmin=" << nmin << ", max_depth=" << max_depth << std::endl;
        recursion_count--;
        return plainmc(f, a, b, N);
    }
    
    // Check for tiny domain size - prevent infinite subdivision
    double domain_size = 0;
    for (int i = 0; i < dim; i++) {
        domain_size += (b[i] - a[i]);
    }
    if (domain_size < 1e-10) {
        std::cout << "Domain too small, terminating recursion" << std::endl;
        recursion_count--;
        return plainmc(f, a, b, N);
    }

    int N_remaining = N - nmin;
    if (N_remaining < 2 * nmin) {
        // Not enough points for stratification
        std::cout << "Not enough points for stratification, using plainmc" << std::endl;
        recursion_count--;
        return plainmc(f, a, b, N);
    }

    // Sample nmin points and evaluate f
    pp::matrix X = sampler(nmin, a, b); // shape: nmin × dim
    pp::vector fxs(nmin);
    for (int i = 0; i < nmin; ++i) fxs[i] = f(X[i]);

    // Find dimension with largest sub-variance
    int best_dim = 0;
    double max_subvar = -1.0;
    
    // Pre-allocate vectors for efficiency
    std::vector<int> left_indices, right_indices;
    left_indices.reserve(nmin);
    right_indices.reserve(nmin);
    
    for (int d = 0; d < dim; ++d) {
        double mid = 0.5 * (a[d] + b[d]);
        
        // Clear previous indices
        left_indices.clear();
        right_indices.clear();
        
        // Partition indices based on dimension d
        for (int i = 0; i < nmin; ++i) {
            if (X(i, d) <= mid)
                left_indices.push_back(i);
            else
                right_indices.push_back(i);
        }
        
        // Skip if any partition is empty
        if (left_indices.empty() || right_indices.empty()) continue;
        
        // Calculate means for each partition
        double sum_left = 0.0, sum_right = 0.0;
        for (int idx : left_indices) sum_left += fxs[idx];
        for (int idx : right_indices) sum_right += fxs[idx];
        
        double muL = sum_left / left_indices.size();
        double muR = sum_right / right_indices.size();
        double subvar = 0.5 * (muL - muR) * (muL - muR);
        
        if (subvar > max_subvar) {
            max_subvar = subvar;
            best_dim = d;
        }
    }
    
    // If we couldn't find a good dimension to split, fall back to plain MC
    if (max_subvar < 0) {
        std::cout << "No good dimension to split, using plainmc" << std::endl;
        recursion_count--;
        return plainmc(f, a, b, N);
    }

    // Split domain along best_dim
    double mid = 0.5 * (a[best_dim] + b[best_dim]);
    pp::vector left_b = b.copy();
    left_b[best_dim] = mid;
    pp::vector right_a = a.copy();
    right_a[best_dim] = mid;

    // Calculate variances for each partition
    left_indices.clear();
    right_indices.clear();
    for (int i = 0; i < nmin; ++i) {
        if (X(i, best_dim) <= mid)
            left_indices.push_back(i);
        else
            right_indices.push_back(i);
    }
    
    // Calculate variances
    double sum_left = 0.0, sum_right = 0.0;
    double sum_sq_left = 0.0, sum_sq_right = 0.0;
    
    for (int idx : left_indices) {
        sum_left += fxs[idx];
        sum_sq_left += fxs[idx] * fxs[idx];
    }
    
    for (int idx : right_indices) {
        sum_right += fxs[idx];
        sum_sq_right += fxs[idx] * fxs[idx];
    }
    
    double n_left = left_indices.size();
    double n_right = right_indices.size();
    
    // Bessel's correction for sample variance
    double vL = n_left > 1 ? (sum_sq_left - (sum_left * sum_left) / n_left) / (n_left - 1) : 0;
    double vR = n_right > 1 ? (sum_sq_right - (sum_right * sum_right) / n_right) / (n_right - 1) : 0;
    
    // Protect against zero variances
    if (vL < 1e-10) vL = 1e-10;
    if (vR < 1e-10) vR = 1e-10;
    
    // Allocate samples based on variance
    int N_left = static_cast<int>(N_remaining * (vL / (vL + vR)));
    int N_right = N_remaining - N_left;
    
    // Ensure minimums
    if (N_left < nmin) {
        N_left = nmin;
        N_right = N_remaining - nmin;
    }
    if (N_right < nmin) {
        N_right = nmin;
        N_left = N_remaining - nmin;
    }
    
    // Check if allocation is possible
    if (N_left < nmin || N_right < nmin) {
        std::cout << "Cannot allocate enough points, using plainmc" << std::endl;
        recursion_count--;
        return plainmc(f, a, b, N);
    }
    
    // Recursive calls with decreased max_depth
    pp::vector left_result = stratifiedmc(f, a, left_b, N_left, nmin, max_depth - 1);
    pp::vector right_result = stratifiedmc(f, right_a, b, N_right, nmin, max_depth - 1);
    
    // Combine results
    double I_total = left_result[0] + right_result[0];
    double err_total = std::sqrt(left_result[1] * left_result[1] + right_result[1] * right_result[1]);
    
    recursion_count--;
    return {I_total, err_total};
}




double test_function(pp::vector x){
    return 1/(1-std::cos(x[0])*std::cos(x[1])*std::cos(x[2]))*1/(M_PI*M_PI*M_PI);
}



double test_function_error(pp::vector x){
    return x[0]*x[1]*x[2]/10.0;
}

double unit_circle(pp::vector x){
    return 1/(1-std::cos(x[0])*std::cos(x[1]))*1/(M_PI*M_PI);
}

double unit_circle_area = M_PI;



// int main(){

//     pp::vector result_list(0);
//     pp::vector deviation_list(0);
//     std::vector N_s = {10, 100, 1000, 10000, 100000, 1000000};

//     pp::vector a = {0, 5, 10};
//     pp::vector b = {5, 10, 15};
//     std::cout << "Calculating deviation for simpler function of type x[0]*x[1]*x[2]/10.0, and creating graph of deviation" << std::endl;
//     for (double N : N_s){

//         pp::vector result = plainmc(test_function_error, a, b, N);
//         result_list.append(result[0]);
//         deviation_list.append(result[1]);
//     }
//     pp::vector::write(result_list, "result_list.txt");
//     pp::vector::write(deviation_list, "deviation_list.txt");
//     // write N_s to file:
//     std::ofstream N_s_file("N_s.txt");
//     for (double N : N_s){
//         N_s_file << N << std::endl;
//     }
//     N_s_file.close();

    
//     std::cout << "Test of plainmc vs quasimc for the difficult singular integral." << std::endl;

//     int N = 1000000;
//     a = {0, 0, 0};
//     b = {M_PI, M_PI, M_PI};
//     pp::vector result = plainmc(test_function, a, b, N);
//     std::cout << "Integral of 1/(1-cos(x)*cos(y)*cos(z)) over [0,pi]^3" << std::endl;
//     std::cout << "Result of integration: " << result[0] << std::endl;
//     std::cout << "Standard deviation: " << result[1] << std::endl;
//     std::cout << "Number of evaluations: " << N << std::endl;
    
//     std::cout << "Test of quasimc" << std::endl;
//     result = quasimc(test_function, a, b, N);
//     std::cout << "Result of integration: " << result[0] << std::endl;
//     std::cout << "Standard deviation: " << result[1] << std::endl;
//     std::cout << "Number of evaluations: " << N << std::endl;


//     std::cout << "Test of stratifiedmc" << std::endl;
//     result = stratifiedmc(test_function_error, a, b, N, 10000);
//     std::cout << "Result of integration: " << result[0] << std::endl;
//     std::cout << "Standard deviation: " << result[1] << std::endl;
//     std::cout << "Number of evaluations: " << N << std::endl;



//     return 0;
// }
int main(){
    pp::vector result_list(0);
    pp::vector deviation_list(0);
    std::vector N_s = {10, 100, 1000, 10000, 100000, 1000000};

    pp::vector a = {0, 5, 10};
    pp::vector b = {5, 10, 15};
    std::cout << "Calculating deviation for simpler function of type x[0]*x[1]*x[2]/10.0, and creating graph of deviation" << std::endl;
    for (double N : N_s){
        pp::vector result = plainmc(test_function_error, a, b, N);
        result_list.append(result[0]);
        deviation_list.append(result[1]);
    }
    pp::vector::write(result_list, "result_list.txt");
    pp::vector::write(deviation_list, "deviation_list.txt");
    // write N_s to file:
    std::ofstream N_s_file("N_s.txt");
    for (double N : N_s){
        N_s_file << N << std::endl;
    }
    N_s_file.close();

    std::cout << "Test of plainmc vs quasimc for the difficult singular integral." << std::endl;

    int N = 1000000;
    a = {0, 0, 0};
    b = {M_PI, M_PI, M_PI};
    pp::vector result = plainmc(test_function, a, b, N);
    std::cout << "Integral of 1/(1-cos(x)*cos(y)*cos(z)) over [0,pi]^3" << std::endl;
    std::cout << "Result of integration: " << result[0] << std::endl;
    std::cout << "Standard deviation: " << result[1] << std::endl;
    std::cout << "Number of evaluations: " << N << std::endl;
    
    std::cout << "Test of quasimc" << std::endl;
    result = quasimc(test_function, a, b, N);
    std::cout << "Result of integration: " << result[0] << std::endl;
    std::cout << "Standard deviation: " << result[1] << std::endl;
    std::cout << "Number of evaluations: " << N << std::endl;

    // // Test stratified MC with smaller N and larger nmin for debugging
    // std::cout << "Test of stratifiedmc with smaller parameters for debugging" << std::endl;
    // // Start with smaller N and larger nmin to debug
    // int test_N = 10000; // Much smaller N
    // int test_nmin = 1000; // Keep nmin relatively large
    // result = stratifiedmc(test_function_error, a, b, test_N, test_nmin, 2); // Limit max depth to 5
    // std::cout << "Result of integration: " << result[0] << std::endl;
    // std::cout << "Standard deviation: " << result[1] << std::endl;
    // std::cout << "Number of evaluations: " << test_N << std::endl;

    // // If that works, try with larger N
    // std::cout << "Now testing with full parameters if first test succeeded" << std::endl;
    // result = stratifiedmc(test_function_error, a, b, N, 1000, 10);
    // std::cout << "Result of integration: " << result[0] << std::endl;
    // std::cout << "Standard deviation: " << result[1] << std::endl;
    // std::cout << "Number of evaluations: " << N << std::endl;

    return 0;
}








