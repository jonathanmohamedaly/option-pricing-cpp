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
        std::cerr << "Error\n";
        return 1;
    }

    out << "nSimulations;Gamma_MC;Gamma_BS;Relative_Error\n";

    double h = 0.01 * S0;

    for (auto n : nSimulations) {
        MonteCarloPricer mc(r, sigma);

        double gamma_mc = mc.gamma(call, S0, T, n, h);
        double gamma_bs = BlackScholes::gamma(S0, K, T, r, sigma);

        double error = std::abs(gamma_mc - gamma_bs) / std::abs(gamma_bs);

        out << n << ";"
            << gamma_mc << ";"
            << gamma_bs << ";"
            << error << "\n";
    }

    out.close();

    std::ofstream out2("convergence_gamma_antithetic.txt");
    if (!out2) {
        std::cerr << "Error\n";
        return 1;
    }

    out2 << "nSimulations;Gamma_Antithetic;Gamma_BS;Relative_Error\n";

    for (auto n : nSimulations) {
        MonteCarloPricer mc(r, sigma);

        double gamma_antithetic = mc.gamma_antithetic(call, S0, T, n, h);
        double gamma_bs = BlackScholes::gamma(S0, K, T, r, sigma);

        double error = std::abs(gamma_antithetic - gamma_bs) / std::abs(gamma_bs);

        out2 << n << ";"
            << gamma_antithetic << ";"
            << gamma_bs << ";"
            << error << "\n";
    }

    out2.close();

    int nSim = 20000;
    int nRuns = 500;

    std::ofstream out3("gamma_histogram.txt");
    out3 << "Run;Gamma_MC;Gamma_Antithetic\n";

    for (int i = 0; i < nRuns; ++i) {
        MonteCarloPricer mc(r, sigma);

        double gamma_mc = mc.gamma(call, S0, T, nSim, h);
        double gamma_ant = mc.gamma_antithetic(call, S0, T, nSim, h);

        out3 << i << ";"
            << gamma_mc << ";"
            << gamma_ant << "\n";
    }

    out3.close();

    nSim = 50000;
    double hFactor = 0.01;

    std::ofstream out4("gamma_vs_spot.txt");
    out4 << "S0;Gamma_MC;Gamma_Antithetic;Gamma_BS\n";

    for (double S0 = 50.0; S0 <= 150.0; S0 += 2.0) {

        double h = hFactor * S0;
        MonteCarloPricer mc(r, sigma);

        double gamma_mc  = mc.gamma(call, S0, T, nSim, h);
        double gamma_ant = mc.gamma_antithetic(call, S0, T, nSim, h);
        double gamma_bs  = BlackScholes::gamma(S0, K, T, r, sigma);

        out4 << S0 << ";"
            << gamma_mc << ";"
            << gamma_ant << ";"
            << gamma_bs << "\n";
    }

    out4.close();

    nSim = 15000;
    nRuns = 300;
    h = 0.01 * S0;

    std::ofstream out5("gamma_boxplot.txt");
    out5 << "Run;Gamma_MC;Gamma_Antithetic\n";

    for (int i = 0; i < nRuns; ++i) {
        MonteCarloPricer mc(r, sigma);

        double gamma_mc  = mc.gamma(call, S0, T, nSim, h);
        double gamma_ant = mc.gamma_antithetic(call, S0, T, nSim, h);

        out5 << i << ";"
            << gamma_mc << ";"
            << gamma_ant << "\n";
    }

    out5.close();

    return 0;
}
