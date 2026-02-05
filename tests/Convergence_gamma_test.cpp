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
    
std::ofstream out("convergence_gamma_mc.txt");
    if (!out) {
        std::cerr << "Erreur : impossible d'ouvrir le fichier\n";
        return 1;
}

out << "nSimulations;Gamma_MC;Gamma_MC_Antithetic;Gamma_BS;Relative_Error\n";

double h = 0.01 * S0;

for (auto n : nSimulations) {
    MonteCarloPricer mc(r, sigma);

    double gamma_mc = mc.gamma(call, S0, T, n, h);
    double gamma_mc_antithetic = mc.gamma_antithetic(call, S0, T, n, h);
    double gamma_bs = BlackScholes::gamma(S0, K, T, r, sigma);

    double error = std::abs(gamma_mc - gamma_bs) / std::abs(gamma_bs);

    out << n << ";"
        << gamma_mc << ";"
        << gamma_mc_antithetic << ";"
        << gamma_bs << ";"
        << error << "\n";
}

out.close();


    return 0;
}
