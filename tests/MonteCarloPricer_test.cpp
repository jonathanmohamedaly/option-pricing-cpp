#include <iostream>
#include "../include/MonteCarloPricer.hpp"
#include "../include/EuropeanOption.hpp"
#include "../include/BlackScholes.hpp"

int main() {
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double sigma = 0.2;

    EuropeanCall call(K);
    EuropeanPut put(K);

    MonteCarloPricer pricer(r, sigma);

    double mcCall = pricer.priceEuropean(S0, T, call, 1000000);
    double bsCall = BlackScholes::call(S0, K, T, r, sigma);

    std::cout << "Monte Carlo Call: " << mcCall << std::endl;
    std::cout << "Black-Scholes Call: " << bsCall << std::endl;
}
