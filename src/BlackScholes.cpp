#include "../include/BlackScholes.hpp"
#include <cmath>

double BlackScholes::norm_cdf(double x){

    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

double BlackScholes::d1(double S_0, double K, double T, double r, double sigma){
    
    return (std::log(S_0 / K) + (r + sigma * sigma / 2) * T) / ( sigma * std::sqrt(T));
}

double BlackScholes::d2(double S_0, double K, double T, double r, double sigma){

    return d1(S_0, K, T, r, sigma) - sigma * std::sqrt(T);
}

double BlackScholes::call(double S_0, double K, double T, double r, double sigma){

    double d1_ = d1(S_0, K, T, r, sigma);
    double d2_ = d2(S_0, K, T, r, sigma);

    return S_0 * norm_cdf(d1_) - K * std::exp(- r * T) * norm_cdf(d2_);
}

double BlackScholes::put(double S_0, double K, double T, double r, double sigma){
    
    double d1_ = d1(S_0, K, T, r, sigma);
    double d2_ = d2(S_0, K, T, r, sigma);

    return -S_0 * norm_cdf(-d1_) + K * std::exp(- r * T) * norm_cdf(-d2_);
}


