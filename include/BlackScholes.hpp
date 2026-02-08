#pragma once 
#include <cmath>

enum class OptionType {Call, Put};
// Compute the Values of Put and Call option using the Black-Scholes equation
class BlackScholes {

    public:

    //Value of d1 in Black-Scholes formula
    static double d1(double S_0, double K, double T, double r, double sigma);

    //Value of d2 in Black-Scholes formula
    static double d2(double S_0, double K, double T, double r, double sigma);

    //Price of an European Call option
    static double call(double S_0, double K, double T, double r, double sigma);

    //Price of an European Put option
    static double put(double S_0, double K, double T, double r, double sigma);

    //Price of an Asian Geometric Call option
    static double geometricAsianCall(double S_0, double K, double r, double sigma, double T, int nSteps);

    //Price of an Asian Geometric Put option
    static double geometricAsianPut(double S_0, double K, double r, double sigma, double T, int nSteps);

    //Greeks

    static double delta(OptionType Option, double S_0, double K, double T, double r, double sigma);

    static double gamma(double S_0, double K, double T, double r, double sigma);

    static double vega(double S_0, double K, double T, double r, double sigma);

    private:

    //CDF of normal std distribution
    static double norm_cdf(double x);
};