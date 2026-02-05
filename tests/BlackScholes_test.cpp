#include "../include/BlackScholes.hpp"
#include <iostream>

int main(){

    double S_0 = 100; //Value of the asset
    double K = 100; //Strike
    double T = 1.0; //Maturity
    double r = 0.05; //Risk-free rate
    double sigma = 0.2; //Volatility

    double call_price = BlackScholes::call(S_0, K, T, r, sigma);
    double put_price = BlackScholes::put(S_0, K, T, r, sigma);

    std::cout << "Call price: " << call_price << std::endl;
    std::cout << "Put price: " << put_price << std::endl;

    double delta_call = BlackScholes::delta(OptionType::Call, S_0, K, T, r, sigma);
    double delta_put  = BlackScholes::delta(OptionType::Put,S_0, K, T, r, sigma);
    double gamma = BlackScholes::gamma(S_0, K, T, r, sigma);
    double vega = BlackScholes::vega(S_0, K, T, r, sigma);

    std::cout << "\n=== Greeks (Black-Scholes) ===\n";
    std::cout << "Delta Call : " << delta_call << std::endl;
    std::cout << "Delta Put  : " << delta_put  << std::endl;
    std::cout << "Gamma      : " << gamma      << std::endl;
    std::cout << "Vega       : " << vega       << std::endl;

    return 0;
}