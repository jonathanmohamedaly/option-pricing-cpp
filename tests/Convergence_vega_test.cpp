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

std::ofstream out("convergence_vega_mc.txt");
if (!out) {
    std::cerr << "Erreur : impossible d'ouvrir le fichier\n";
    return 1;
}

out << "nSimulations;Vega_MC;Vega_MC_Antithetic;Vega_BS;Relative_Error\n";

double h = 0.01;

for (auto n : nSimulations) {
    MonteCarloPricer mc(r, sigma);

    double vega_mc = mc.vega(call, S0, T, n, h);
    double vega_mc_antithetic = mc.vega_antithetic(call, S0, T, n, h);
    double vega_bs = BlackScholes::vega(S0, K, T, r, sigma);

    double error = std::abs(vega_mc - vega_bs) / std::abs(vega_bs);

    out << n << ";"
        << vega_mc << ";"
        << vega_mc_antithetic << ";"
        << vega_bs << ";"
        << error << "\n";
}

out.close();



    return 0;
}
