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
    
    std::ofstream out("convergence_delta_mc.txt");
if (!out) {
    std::cerr << "Erreur : impossible d'ouvrir le fichier\n";
    return 1;
}

out << "nSimulations;Delta_MC;Delta_MC_Antithetic;Delta_BS;Relative_Error\n";

double h = 0.01 * S0;

for (auto n : nSimulations) {
    MonteCarloPricer mc(r, sigma);

    double delta_mc = mc.delta(call, S0, T, n, h);
    double delta_mc_antithetic = mc.delta_antithetic(call, S0, T, n, h);
    double delta_bs = BlackScholes::delta(OptionType::Call, S0, K, T, r, sigma);

    double error = std::abs(delta_mc - delta_bs) / std::abs(delta_bs);

    out << n << ";"
        << delta_mc << ";"
        << delta_mc_antithetic << ";"
        << delta_bs << ";"
        << error << "\n";
}

out.close();

    return 0;
}
