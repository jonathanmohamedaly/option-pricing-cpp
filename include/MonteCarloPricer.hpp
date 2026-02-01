#pragma once 

#include "Payoff.hpp"
#include "RandomGenerator.hpp"


//Estimate the value of an option using Monte Carlo Method
class MonteCarloPricer {

    public:

        MonteCarloPricer(double r, double sigma);

        double priceEuropean(double S_0, double T, const Payoff& payoff, int nSimulations) const;

    private:

        double r_;
        double sigma_;


};