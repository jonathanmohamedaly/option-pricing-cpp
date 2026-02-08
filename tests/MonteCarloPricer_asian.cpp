#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

#include "../include/MonteCarloPricer.hpp"
#include "../include/AsianOption.hpp"

int main() {

    double S0 = 100.0;
    double K  = 100.0;
    double r  = 0.05;
    double sigma = 0.2;
    double T  = 1.0;

    int nSteps = 50;

    AsianCall arithCall(K);
    AsianGeometricCall  geomCall(K);

    int nMin = 1'000;
    int nMax = 200'000;
    int nPoints = 100;

    std::vector<int> nSimulations;
    for (int i = 0; i < nPoints; ++i) {
        int n = nMin + i * (nMax - nMin) / (nPoints - 1);
        nSimulations.push_back(n);
    }

    MonteCarloPricer mc(r, sigma);

    // MC
    int nRef = 1'000'000;
    double referencePrice = mc.priceAsian(S0, T, arithCall, nRef, nSteps);

    // MC antithetic
    std::ofstream out2("convergence_mc_antithetic_asian.txt");
    out2 << "nSimulations;Price;Relative_Error\n";

    for (int n : nSimulations) {
        double price = mc.priceAsianAntithetic(S0, T, arithCall, n, nSteps);
        double error = std::abs(price - referencePrice) / referencePrice;

        out2 << n << ";" << price << ";" << error << "\n";
    }
    out2.close();

    // Control variate 
    std::ofstream out3("convergence_mc_control_variate_asian.txt");
    out3 << "nSimulations;Price;Relative_Error\n";

    for (int n : nSimulations) {
        double price = mc.priceAsianControlVariate(S0, K, OptionType::Call, T, arithCall, n, nSteps);

        double error = std::abs(price - referencePrice) / referencePrice;

        out3 << n << ";" << price << ";" << error << "\n";
    }
    out3.close();

    // Control variate géométrique 
    std::ofstream out4("convergence_mc_control_variate_geometric.txt");
    out4 << "nSimulations;Price;Relative_Error\n";

    for (int n : nSimulations) {
        double price = mc.priceAsianGeometricControlVariate(S0, K, OptionType::Call, T, arithCall, geomCall, n, nSteps);

        double error = std::abs(price - referencePrice) / referencePrice;

        out4 << n << ";" << price << ";" << error << "\n";
    }
    out4.close();


    return 0;
}
