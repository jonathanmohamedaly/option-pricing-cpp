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
    
    std::ofstream out("convergence_vega_mc.txt");
    if (!out) {
        std::cerr << "Error\n";
        return 1;
    }

    out << "nSimulations;Vega_MC;Vega_BS;Relative_Error\n";

    double h = 0.075 * S0;

    for (auto n : nSimulations) {
        MonteCarloPricer mc(r, sigma);

        double vega_mc = mc.vega(call, S0, T, n, h);
        double vega_bs = BlackScholes::vega(S0, K, T, r, sigma);

        double error = std::abs(vega_mc - vega_bs) / std::abs(vega_bs);

        out << n << ";"
            << vega_mc << ";"
            << vega_bs << ";"
            << error << "\n";
    }

    out.close();

    std::ofstream out2("convergence_vega_antithetic.txt");
    if (!out2) {
        std::cerr << "Error\n";
        return 1;
    }

    out2 << "nSimulations;Vega_Antithetic;Vega_BS;Relative_Error\n";

    for (auto n : nSimulations) {
        MonteCarloPricer mc(r, sigma);

        double vega_antithetic = mc.vega_antithetic(call, S0, T, n, h);
        double vega_bs = BlackScholes::vega(S0, K, T, r, sigma);

        double error = std::abs(vega_antithetic - vega_bs) / std::abs(vega_bs);

        out2 << n << ";"
            << vega_antithetic << ";"
            << vega_bs << ";"
            << error << "\n";
    }

    out2.close();

    int nSim = 20000;
    int nRuns = 500;

    std::ofstream out3("vega_histogram.txt");
    out3 << "Run;Vega_MC;Vega_Antithetic\n";

    for (int i = 0; i < nRuns; ++i) {
        MonteCarloPricer mc(r, sigma);

        double vega_mc = mc.vega(call, S0, T, nSim, h);
        double vega_ant = mc.vega_antithetic(call, S0, T, nSim, h);

        out3 << i << ";"
            << vega_mc << ";"
            << vega_ant << "\n";
    }

    out3.close();

    nSim = 50000;
    double hFactor = 0.075;

    std::ofstream out4("vega_vs_spot.txt");
    out4 << "S0;Vega_MC;Vega_Antithetic;Vega_BS\n";

    for (double S0 = 50.0; S0 <= 150.0; S0 += 2.0) {

        double h = hFactor * S0;
        MonteCarloPricer mc(r, sigma);

        double vega_mc  = mc.vega(call, S0, T, nSim, h);
        double vega_ant = mc.vega_antithetic(call, S0, T, nSim, h);
        double vega_bs  = BlackScholes::vega(S0, K, T, r, sigma);

        out4 << S0 << ";"
            << vega_mc << ";"
            << vega_ant << ";"
            << vega_bs << "\n";
    }

    out4.close();

    nSim = 15000;
    nRuns = 300;
    h = 0.075 * S0;

    std::ofstream out5("vega_boxplot.txt");
    out5 << "Run;Vega_MC;Vega_Antithetic\n";

    for (int i = 0; i < nRuns; ++i) {
        MonteCarloPricer mc(r, sigma);

        double vega_mc  = mc.vega(call, S0, T, nSim, h);
        double vega_ant = mc.vega_antithetic(call, S0, T, nSim, h);

        out5 << i << ";"
            << vega_mc << ";"
            << vega_ant << "\n";
    }

    out5.close();

    return 0;
}
