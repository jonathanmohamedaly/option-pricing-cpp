#pragma once
#include <vector>

class Payoff {

    public:

    virtual ~Payoff() = default;

    virtual double operator()(double ST) const = 0;

    virtual double operator()(const std::vector<double>& path) const {
        return (*this)(path.back());
    }
};