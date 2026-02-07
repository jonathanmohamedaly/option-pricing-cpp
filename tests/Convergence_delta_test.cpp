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
        std::cerr << "Error\n";
        return 1;
    }

    out << "nSimulations;Delta_MC;Delta_BS;Relative_Error\n";

    double h = 0.01 * S0;

    for (auto n : nSimulations) {
        MonteCarloPricer mc(r, sigma);

        double delta_mc = mc.delta(call, S0, T, n, h);
        double delta_bs = BlackScholes::delta(OptionType::Call, S0, K, T, r, sigma);

        double error = std::abs(delta_mc - delta_bs) / std::abs(delta_bs);

        out << n << ";"
            << delta_mc << ";"
            << delta_bs << ";"
            << error << "\n";
    }

    out.close();

    std::ofstream out2("convergence_delta_antithetic.txt");
    if (!out2) {
        std::cerr << "Error\n";
        return 1;
    }

    out2 << "nSimulations;Delta_Antithetic;Delta_BS;Relative_Error\n";

    for (auto n : nSimulations) {
        MonteCarloPricer mc(r, sigma);

        double delta_antithetic = mc.delta_antithetic(call, S0, T, n, h);
        double delta_bs = BlackScholes::delta(OptionType::Call, S0, K, T, r, sigma);

        double error = std::abs(delta_antithetic - delta_bs) / std::abs(delta_bs);

        out2 << n << ";"
            << delta_antithetic << ";"
            << delta_bs << ";"
            << error << "\n";
    }

    out2.close();

    int nSim = 20000;
    int nRuns = 500;

    std::ofstream out3("delta_histogram.txt");
    out3 << "Run;Delta_MC;Delta_Antithetic\n";

    for (int i = 0; i < nRuns; ++i) {
        MonteCarloPricer mc(r, sigma);

        double delta_mc = mc.delta(call, S0, T, nSim, h);
        double delta_ant = mc.delta_antithetic(call, S0, T, nSim, h);

        out3 << i << ";"
            << delta_mc << ";"
            << delta_ant << "\n";
    }

    out3.close();

    nSim = 50000;
    double hFactor = 0.01;

    std::ofstream out4("delta_vs_spot.txt");
    out4 << "S0;Delta_MC;Delta_Antithetic;Delta_BS\n";

    for (double S0 = 50.0; S0 <= 150.0; S0 += 2.0) {

        double h = hFactor * S0;
        MonteCarloPricer mc(r, sigma);

        double delta_mc  = mc.delta(call, S0, T, nSim, h);
        double delta_ant = mc.delta_antithetic(call, S0, T, nSim, h);
        double delta_bs  = BlackScholes::delta(OptionType::Call, S0, K, T, r, sigma);

        out4 << S0 << ";"
            << delta_mc << ";"
            << delta_ant << ";"
            << delta_bs << "\n";
    }

    out4.close();

    nSim = 15000;
    nRuns = 300;
    h = 0.01 * S0;

    std::ofstream out5("delta_boxplot.txt");
    out5 << "Run;Delta_MC;Delta_Antithetic\n";

    for (int i = 0; i < nRuns; ++i) {
        MonteCarloPricer mc(r, sigma);

        double delta_mc  = mc.delta(call, S0, T, nSim, h);
        double delta_ant = mc.delta_antithetic(call, S0, T, nSim, h);

        out5 << i << ";"
            << delta_mc << ";"
            << delta_ant << "\n";
    }

    out5.close();

    return 0;
}
