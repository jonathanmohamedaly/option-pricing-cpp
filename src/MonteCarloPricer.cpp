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

//Useless for European Options
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

        // Discount directement ici
        double X1 = std::exp(-r_*T) * payoff(S_T1);
        double X2 = std::exp(-r_*T) * payoff(S_T2);

        double Y1, Y2;
        if (type == OptionType::Call) {
            Y1 = std::exp(-r_*T) * std::max(S_T1 - K, 0.0);
            Y2 = std::exp(-r_*T) * std::max(S_T2 - K, 0.0);
        } else {
            Y1 = std::exp(-r_*T) * std::max(K - S_T1, 0.0);
            Y2 = std::exp(-r_*T) * std::max(K - S_T2, 0.0);
        }

        sum_X  += X1 + X2;
        sum_Y  += Y1 + Y2;
        sum_XY += X1*Y1 + X2*Y2;
        sum_YY += Y1*Y1 + Y2*Y2;
    }

    int totalSamples = 2 * nSimulations;

    double mean_X = sum_X / totalSamples;
    double mean_Y = sum_Y / totalSamples;

    double cov_XY = sum_XY / totalSamples - mean_X * mean_Y;
    double var_Y  = sum_YY / totalSamples - mean_Y * mean_Y;

    double beta = cov_XY / var_Y;

    double price = mean_X + beta * (bs_price - mean_Y);

    return price;
}



double MonteCarloPricer::delta(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const {

    double payoff_plus  = 0.0;
    double payoff_minus = 0.0;

    for (int i = 0; i < nSimulations; ++i) {
        double Z = rng_.gaussian();

        double drift = (r_ - 0.5 * sigma_ * sigma_) * T;
        double vol   = sigma_ * std::sqrt(T);

        double S_plus  = (S_0 + h) * std::exp(drift + vol * Z);
        double S_minus = (S_0 - h) * std::exp(drift + vol * Z);

        payoff_plus  += payoff(S_plus);
        payoff_minus += payoff(S_minus);
    }

    double price_plus  = std::exp(-r_ * T) * payoff_plus  / nSimulations;
    double price_minus = std::exp(-r_ * T) * payoff_minus / nSimulations;

    return (price_plus - price_minus) / (2.0 * h);
}

double MonteCarloPricer::delta_antithetic(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const {
    
    if (nSimulations % 2 != 0) nSimulations++;

    double payoff_plus  = 0.0;
    double payoff_minus = 0.0;

    for (int i = 0; i < nSimulations / 2; ++i) {
        double Z = rng_.gaussian();

        double drift = (r_ - 0.5 * sigma_ * sigma_) * T;
        double vol   = sigma_ * std::sqrt(T);

        double S_p1 = (S_0 + h) * std::exp(drift + vol * Z);
        double S_p2 = (S_0 + h) * std::exp(drift - vol * Z);

        double S_m1 = (S_0 - h) * std::exp(drift + vol * Z);
        double S_m2 = (S_0 - h) * std::exp(drift - vol * Z);

        payoff_plus  += payoff(S_p1) + payoff(S_p2);
        payoff_minus += payoff(S_m1) + payoff(S_m2);
    }

    payoff_plus  /= nSimulations;
    payoff_minus /= nSimulations;

    double price_plus  = std::exp(-r_ * T) * payoff_plus;
    double price_minus = std::exp(-r_ * T) * payoff_minus;

    return (price_plus - price_minus) / (2.0 * h);
}

double MonteCarloPricer::gamma(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const {
    
    double payoff_plus  = 0.0;
    double payoff_mid   = 0.0;
    double payoff_minus = 0.0;

    for (int i = 0; i < nSimulations; ++i) {
        double Z = rng_.gaussian();

        double drift = (r_ - 0.5 * sigma_ * sigma_) * T;
        double vol   = sigma_ * std::sqrt(T);

        double ST_plus  = (S_0 + h) * std::exp(drift + vol * Z);
        double ST_mid   = S_0 * std::exp(drift + vol * Z);
        double ST_minus = (S_0 - h) * std::exp(drift + vol * Z);

        payoff_plus  += payoff(ST_plus);
        payoff_mid   += payoff(ST_mid);
        payoff_minus += payoff(ST_minus);
    }

    payoff_plus  /= nSimulations;
    payoff_mid   /= nSimulations;
    payoff_minus /= nSimulations;

    double disc = std::exp(-r_ * T);

    return disc * (payoff_plus - 2.0 * payoff_mid + payoff_minus) / (h * h);
}

double MonteCarloPricer::gamma_antithetic(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const {

    if (nSimulations % 2 != 0) nSimulations++;

    double payoff_plus  = 0.0;
    double payoff_mid   = 0.0;
    double payoff_minus = 0.0;

    for (int i = 0; i < nSimulations / 2; ++i) {
        double Z = rng_.gaussian();

        double drift = (r_ - 0.5 * sigma_ * sigma_) * T;
        double vol   = sigma_ * std::sqrt(T);

        double ST_p1 = (S_0 + h) * std::exp(drift + vol * Z);
        double ST_p2 = (S_0 + h) * std::exp(drift - vol * Z);

        double ST_m1 = (S_0 - h) * std::exp(drift + vol * Z);
        double ST_m2 = (S_0 - h) * std::exp(drift - vol * Z);

        double ST_0_1 = S_0 * std::exp(drift + vol * Z);
        double ST_0_2 = S_0 * std::exp(drift - vol * Z);

        payoff_plus  += payoff(ST_p1) + payoff(ST_p2);
        payoff_minus += payoff(ST_m1) + payoff(ST_m2);
        payoff_mid   += payoff(ST_0_1) + payoff(ST_0_2);
    }

    payoff_plus  /= nSimulations;
    payoff_mid   /= nSimulations;
    payoff_minus /= nSimulations;

    double disc = std::exp(-r_ * T);

    return disc * (payoff_plus - 2.0 * payoff_mid + payoff_minus) / (h * h);
}

double MonteCarloPricer::vega(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const {
    
    double payoff_plus  = 0.0;
    double payoff_minus = 0.0;

    for (int i = 0; i < nSimulations; ++i) {

        double Z = rng_.gaussian();

        double drift_plus  = (r_ - 0.5 * (sigma_ + h) * (sigma_ + h)) * T;
        double drift_minus = (r_ - 0.5 * (sigma_ - h) * (sigma_ - h)) * T;

        double vol_plus  = (sigma_ + h) * std::sqrt(T);
        double vol_minus = (sigma_ - h) * std::sqrt(T);

        double S_plus  = S_0 * std::exp(drift_plus  + vol_plus  * Z);
        double S_minus = S_0 * std::exp(drift_minus + vol_minus * Z);

        payoff_plus  += payoff(S_plus);
        payoff_minus += payoff(S_minus);
    }

    payoff_plus  /= nSimulations;
    payoff_minus /= nSimulations;

    double price_plus  = std::exp(-r_ * T) * payoff_plus;
    double price_minus = std::exp(-r_ * T) * payoff_minus;

    return (price_plus - price_minus) / (2.0 * h);
}

double MonteCarloPricer::vega_antithetic(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const {
    
    if (nSimulations % 2 != 0) nSimulations++;

    double payoff_plus  = 0.0;
    double payoff_minus = 0.0;

    for (int i = 0; i < nSimulations / 2; ++i) {

        double Z = rng_.gaussian();

        double drift_plus  = (r_ - 0.5 * (sigma_ + h) * (sigma_ + h)) * T;
        double drift_minus = (r_ - 0.5 * (sigma_ - h) * (sigma_ - h)) * T;

        double vol_plus  = (sigma_ + h) * std::sqrt(T);
        double vol_minus = (sigma_ - h) * std::sqrt(T);

        double S_p1 = S_0 * std::exp(drift_plus  + vol_plus  * Z);
        double S_p2 = S_0 * std::exp(drift_plus  - vol_plus  * Z);

        double S_m1 = S_0 * std::exp(drift_minus + vol_minus * Z);
        double S_m2 = S_0 * std::exp(drift_minus - vol_minus * Z);

        payoff_plus  += payoff(S_p1) + payoff(S_p2);
        payoff_minus += payoff(S_m1) + payoff(S_m2);
    }

    payoff_plus  /= nSimulations;
    payoff_minus /= nSimulations;

    double price_plus  = std::exp(-r_ * T) * payoff_plus;
    double price_minus = std::exp(-r_ * T) * payoff_minus;

    return (price_plus - price_minus) / (2.0 * h);

}