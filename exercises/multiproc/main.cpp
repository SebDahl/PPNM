#include <iostream>
#include <vector>
#include <thread>
#include <string>
#include <sstream>
#include <cstdlib>

struct Data {
    int a, b;
    double sum = 0.0;
};

void harm(Data* arg) {
    arg->sum = 0.0;
    for (int i = arg->a; i < arg->b; ++i)
        arg->sum += 1.0 / i;
}

int main(int argc, char* argv[]) {
    int nthreads = 1;
    int nterms = static_cast<int>(1e8);

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        size_t pos = arg.find(':');
        if (pos != std::string::npos) {
            std::string key = arg.substr(0, pos);
            std::string value = arg.substr(pos + 1);
            if (key == "-threads") nthreads = std::stoi(value);
            if (key == "-terms") nterms = static_cast<int>(std::stof(value));
        }
    }

    // Prepare data objects
    std::vector<Data> params(nthreads);
    for (int i = 0; i < nthreads; ++i) {
        params[i].a = 1 + nterms / nthreads * i;
        params[i].b = 1 + nterms / nthreads * (i + 1);
    }
    params[nthreads - 1].b = nterms + 1; // Adjust the endpoint

    // Create and start threads
    std::vector<std::thread> threads;
    for (int i = 0; i < nthreads; ++i) {
        threads.emplace_back(harm, &params[i]);
    }

    // Join threads
    for (auto& t : threads) t.join();

    // Sum all results
    double total = 0.0;
    for (const auto& p : params)
        total += p.sum;

    std::cout.precision(15);
    std::cout << "Harmonic sum = " << total << std::endl;

    return 0;
}
