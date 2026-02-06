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

double MonteCarloPricer::priceEuropeanControlVariate(double S_0, double K, OptionType type, double T, const Payoff& payoff, int nSimulations) const{

    double sum_X = 0.0, sum_Y = 0.0;
    double sum_XY = 0.0, sum_YY = 0.0;

    double bs_price =
        (type == OptionType::Call)
        ? BlackScholes::call(S_0, K, T, r_, sigma_)
        : BlackScholes::put (S_0, K, T, r_, sigma_);

    for (int i = 0; i < nSimulations; ++i) {
        double Z = rng_.gaussian();

        double S_T1 = S_0 * std::exp((r_ - 0.5*sigma_*sigma_)*T + sigma_*sqrt(T)*Z);
        double S_T2 = S_0 * std::exp((r_ - 0.5*sigma_*sigma_)*T + sigma_*sqrt(T)*(-Z));

        double X1 = payoff(S_T1);
        double X2 = payoff(S_T2);

        double Y1 = std::max(S_T1 - K, 0.0);
        double Y2 = std::max(S_T2 - K, 0.0);

        sum_X  += X1 + X2;
        sum_Y  += Y1 + Y2;
        sum_XY += X1*Y1 + X2*Y2;
        sum_YY += Y1*Y1 + Y2*Y2;
    }

    double mean_X = sum_X / nSimulations;
    double mean_Y = sum_Y / nSimulations;

    double cov_XY = sum_XY / nSimulations - mean_X * mean_Y;
    double var_Y  = sum_YY / nSimulations - mean_Y * mean_Y;

    double beta = cov_XY / var_Y;

    double price = mean_X + beta * (bs_price - mean_Y);

    return std::exp(-r_ * T) * price;
}



double MonteCarloPricer::delta(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const {
    double sup = priceEuropean(S_0 + h, T, payoff, nSimulations);
    double inf = priceEuropean(S_0 - h, T, payoff, nSimulations);

    return (sup - inf) / (2.0 * h);
}

double MonteCarloPricer::delta_antithetic(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const {
    double sup = priceEuropeanAntithetic(S_0 + h, T, payoff, nSimulations);
    double inf = priceEuropeanAntithetic(S_0 - h, T, payoff, nSimulations);

    return (sup - inf) / (2.0 * h);
}

double MonteCarloPricer::gamma(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const {
    double sup = priceEuropean(S_0 + h, T, payoff, nSimulations);
    double mid = priceEuropean(S_0, T, payoff, nSimulations);
    double inf = priceEuropean(S_0 - h, T, payoff, nSimulations);

    return (sup - 2.0 * mid + inf) / (h * h);
}

double MonteCarloPricer::gamma_antithetic(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const {
    double sup = priceEuropeanAntithetic(S_0 + h, T, payoff, nSimulations);
    double mid = priceEuropeanAntithetic(S_0, T, payoff, nSimulations);
    double inf = priceEuropeanAntithetic(S_0 - h, T, payoff, nSimulations);

    return (sup - 2.0 * mid + inf) / (h * h);
}

double MonteCarloPricer::vega(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const {
    MonteCarloPricer mc_sup(r_, sigma_ + h);
    MonteCarloPricer mc_inf(r_, sigma_ - h);

    double sup  = mc_sup.priceEuropean(S_0, T, payoff, nSimulations);
    double inf = mc_inf.priceEuropean(S_0, T, payoff, nSimulations);

    return (sup - inf) / (2.0 * h);
}

double MonteCarloPricer::vega_antithetic(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const {
    MonteCarloPricer mc_sup(r_, sigma_ + h);
    MonteCarloPricer mc_inf(r_, sigma_ - h);

    double sup  = mc_sup.priceEuropeanAntithetic(S_0, T, payoff, nSimulations);
    double inf = mc_inf.priceEuropeanAntithetic(S_0, T, payoff, nSimulations);

    return (sup - inf) / (2.0 * h);
}