#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

#include "../include/MonteCarloPricer.hpp"
#include "../include/AsianOption.hpp"
#include "../include/BlackScholes.hpp"

int main() {

    double S0 = 100.0;
    double K  = 100.0;
    double r  = 0.05;
    double sigma = 0.2;
    double T  = 1.0;

    int nSteps = 50;

    AsianCall arithCall(K);
    AsianGeometricCall geomCall(K);

    MonteCarloPricer mc(r, sigma);

    // Delta exact géométrique pour référence
    double delta_ref =
        BlackScholes::geometricAsianCall(S0, K, r, sigma, T, nSteps);

    // Convergence
    int nMin = 1'000;
    int nMax = 200'000;
    int nPoints = 100;

    std::vector<int> nSimulations;
    for (int i = 0; i < nPoints; ++i) {
        int n = nMin + i * (nMax - nMin) / (nPoints - 1);
        nSimulations.push_back(n);
    }

    std::ofstream out1("convergence_delta_asian_pathwise.txt");
    out1 << "nSim;Delta_PW;Delta_PW_Ant;Delta_PW_Ant_CV;RelErr_PW;RelErr_PW_Ant;RelErr_PW_Ant_CV\n";

    for (int n : nSimulations) {

        double d_pw      = mc.deltaAsianPathwise(S0, K, OptionType::Call, T, n, nSteps);
        double d_pw_ant  = mc.deltaAsianPathwiseAntithetic(S0, K, OptionType::Call, T, n, nSteps);
        double d_pw_cv   = mc.deltaAsianPathwiseAntitheticCV(S0, K, OptionType::Call, T, n, nSteps);

        out1 << n << ";"
             << d_pw << ";"
             << d_pw_ant << ";"
             << d_pw_cv << ";"
             << std::abs(d_pw      - delta_ref)/delta_ref << ";"
             << std::abs(d_pw_ant  - delta_ref)/delta_ref << ";"
             << std::abs(d_pw_cv   - delta_ref)/delta_ref << "\n";
    }

    out1.close();

    // Histograms
    int nSim = 20'000;
    int nRuns = 400;

    std::ofstream out2("delta_asian_histogram_pathwise.txt");
    out2 << "Run;Delta_PW;Delta_PW_Ant;Delta_PW_Ant_CV\n";

    for (int i = 0; i < nRuns; ++i) {

        double d_pw      = mc.deltaAsianPathwise(S0, K, OptionType::Call, T, nSim, nSteps);
        double d_pw_ant  = mc.deltaAsianPathwiseAntithetic(S0, K, OptionType::Call, T, nSim, nSteps);
        double d_pw_cv   = mc.deltaAsianPathwiseAntitheticCV(S0, K, OptionType::Call, T, nSim, nSteps);

        out2 << i << ";"
             << d_pw << ";"
             << d_pw_ant << ";"
             << d_pw_cv << "\n";
    }

    out2.close();

    // Delta vs Spot
    nSim = 50'000;

    std::ofstream out3("delta_asian_vs_spot_pathwise.txt");
    out3 << "S0;Delta_PW;Delta_PW_Ant;Delta_PW_Ant_CV;Delta_Geo_Exact\n";

    for (double spot = 50.0; spot <= 150.0; spot += 2.0) {

        double d_pw      = mc.deltaAsianPathwise(spot, K, OptionType::Call, T, nSim, nSteps);
        double d_pw_ant  = mc.deltaAsianPathwiseAntithetic(spot, K, OptionType::Call, T, nSim, nSteps);
        double d_pw_cv   = mc.deltaAsianPathwiseAntitheticCV(spot, K, OptionType::Call, T, nSim, nSteps);

        double d_ref     = BlackScholes::geometricAsianCall(spot, K, r, sigma, T, nSteps);

        out3 << spot << ";"
             << d_pw << ";"
             << d_pw_ant << ";"
             << d_pw_cv << ";"
             << d_ref << "\n";
    }

    out3.close();

    return 0;
}