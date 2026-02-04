#include "../include/MonteCarloPricer.hpp"
#include <cmath>

MonteCarloPricer::MonteCarloPricer(double r, double sigma){
    r_ = r;
    sigma_ = sigma;
    rng_ = RandomGenerator();
}

double MonteCarloPricer::priceEuropean(double S_0, double T, const Payoff& payoff, int nSimulations) const{

    double total_payoff = 0.0;

    for(int i = 0; i < nSimulations; i++){

        double Z = rng_.gaussian();
        double S_T = S_0 * std::exp((r_ - 0.5 * sigma_ * sigma_) * T + sigma_ * std::sqrt(T) * Z);
        total_payoff += payoff(S_T);

    }

    return std::exp(-r_ * T) * total_payoff / nSimulations;
}

double MonteCarloPricer::priceEuropeanAntithetic(double S_0, double T, const Payoff& payoff, int nSimulations) const {
    if (nSimulations % 2 != 0) nSimulations++;  

    double total_payoff = 0.0;

    for(int i = 0; i < nSimulations / 2; i++){
        double Z = rng_.gaussian();
        double S_T1 = S_0 * std::exp((r_ - 0.5 * sigma_ * sigma_) * T + sigma_ * std::sqrt(T) * Z);
        double S_T2 = S_0 * std::exp((r_ - 0.5 * sigma_ * sigma_) * T + sigma_ * std::sqrt(T) * (-Z));

        total_payoff += payoff(S_T1) + payoff(S_T2);  
    }

    return std::exp(-r_ * T) * total_payoff / nSimulations;  
}
