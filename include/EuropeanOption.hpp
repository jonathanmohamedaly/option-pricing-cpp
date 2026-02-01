#pragma once
#include "Payoff.hpp"

class EuropeanCall : public Payoff {

    public:
        explicit EuropeanCall(double strike) : K_(strike) {}
        double operator()(double ST) const override {
            return std::max(ST - K_, 0.0);
        }
        
    private:
        double K_;
};


class EuropeanPut : public Payoff {

    public:
        explicit EuropeanPut(double strike) : K_(strike) {}
        double operator()(double ST) const override {
            return std::max(K_ - ST, 0.0);
        }
        
    private:
        double K_;
};