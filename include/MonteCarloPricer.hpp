#pragma once 

#include "Payoff.hpp"
#include "RandomGenerator.hpp"
#include "BlackScholes.hpp"

//Estimate the value of an option using Monte Carlo Method
class MonteCarloPricer {

    public:

        MonteCarloPricer(double r, double sigma);

        double priceEuropean(double S_0, double T, const Payoff& payoff, int nSimulations) const;

        double priceEuropeanAntithetic(double S_0, double T, const Payoff& payoff, int nSimulations) const;

        double priceEuropeanControlVariate(double S_0, double K, OptionType type, double T, const Payoff& payoff, int nSimulations) const;

        double priceAsian(double S_0, double T, const Payoff& payoff, int nSimulations, int nSteps) const;

        double priceAsianAntithetic(double S_0, double T, const Payoff& payoff, int nSimulations, int nSteps) const;
        
        double priceAsianControlVariate(double S_0, double K, OptionType type, double T, const Payoff& asianPayoff, int nSimulations,
                                        int nSteps) const;

        double priceAsianGeometricControlVariate(double S_0, double K, OptionType type, double T, const Payoff& asianArithmeticPayoff,
                                                 const Payoff& asianGeometricPayoff, int nSimulations, int nSteps) const;

        double delta(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const;

        double delta_antithetic(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const;

        double gamma(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const;

        double gamma_antithetic(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const;

        double vega(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const;

        double vega_antithetic(const Payoff& payoff, double S_0, double T, int nSimulations, double h) const;

        double deltaAsianPathwise(double S0, double K, OptionType type, double T, int nSimulations, int nSteps) const ;

        double deltaAsianPathwiseAntithetic(double S0, double K, OptionType type, double T, int nSimulations, int nSteps) const;

        double deltaAsianPathwiseAntitheticCV( double S0, double K, OptionType type, double T, int nSimulations, int nSteps) const;

    private:

        double r_;
        double sigma_;
        mutable RandomGenerator rng_;

};