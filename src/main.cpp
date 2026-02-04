#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "../include/MonteCarloPricer.hpp" 
#include "../include/BlackScholes.hpp"     
#include "../include/EuropeanOption.hpp"     

int main() {
    double S0 = 100.0, K = 100.0, r = 0.05, sigma = 0.2, T = 1.0;
    EuropeanCall call(K);

    std::vector<int> nSimulations;
    nSimulations.reserve(500);

    int nMin = 1000;
    int nMax = 1000000;
    int nPoints = 1000;

    for (int i = 0; i < nPoints; ++i) {
        int n = nMin + i * (nMax - nMin) / (nPoints - 1);
        nSimulations.push_back(n);
    }
    double bsPrice = BlackScholes::call(S0, K, T, r, sigma);

    std::ofstream out("convergence_mc.txt");
    if (!out) {
        std::cerr << "Erreur : impossible d'ouvrir le fichier\n";
        return 1;
    }

    out << "nSimulations;MC_Price;Relative_Error\n";

    for (auto n : nSimulations) {
        MonteCarloPricer mc(r, sigma);
        double mcPrice = mc.priceEuropean(S0, T, call, n);
        double error = std::abs(mcPrice - bsPrice) / bsPrice;

        out << n << ";"
            << mcPrice << ";"
            << error << "\n";
    }

    out.close();

    std::ofstream out2("convergence_mc_antithetic.txt");
    if (!out2) {
        std::cerr << "Erreur : impossible d'ouvrir le fichier\n";
        return 1;
    }

    out2 << "nSimulations;MC_Price_Antithetic;Relative_Error\n";

    for (auto n : nSimulations) {
        MonteCarloPricer mc(r, sigma);
        double mcPrice = mc.priceEuropeanAntithetic(S0, T, call, n);
        double error = std::abs(mcPrice - bsPrice) / bsPrice;

        out2 << n << ";"
            << mcPrice << ";"
            << error << "\n";
    }

    out2.close();

    return 0;
}
