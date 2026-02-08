#pragma once
#include "Payoff.hpp"
#include <vector>
#include <algorithm> 
#include <numeric>

class AsianCall : public Payoff {

    public: 
        explicit AsianCall(double strike) : K_(strike){}

        double operator()(const std::vector<double>& path) const override {

            double sum = std::accumulate(path.begin() + 1, path.end(), 0.0);
            double average = sum / (path.size() - 1);
            return std::max(average - K_, 0.0);
        }

        double operator()(double ST) const override {
            return std::max(ST - K_, 0.0); 
        }

    private:
        double K_;
};

class AsianPut : public Payoff {

    public: 
        explicit AsianPut(double strike) : K_(strike){}

        double operator()(const std::vector<double>& path) const override {

            double sum = std::accumulate(path.begin() + 1, path.end(), 0.0);
            double average = sum / (path.size() - 1);
            return std::max(K_-average, 0.0);
        }

        double operator()(double ST) const override {
            return std::max(K_ - ST, 0.0); 
        }

    private:
        double K_;
};

class AsianGeometricCall : public Payoff {

    public:
        explicit AsianGeometricCall(double strike) : K_(strike) {}

        double operator()(const std::vector<double>& path) const override {
            double logSum = 0.0;
            for (size_t i = 1; i < path.size(); ++i)
                logSum += std::log(path[i]);

            double geomMean = std::exp(logSum / (path.size() - 1));
            return std::max(geomMean - K_, 0.0);
        }

        double operator()(double ST) const override {
            return std::max(ST - K_, 0.0); 
        }

    private:
        double K_;
};

class AsianGeometricPut : public Payoff {

    public:
        explicit AsianGeometricPut(double strike) : K_(strike) {}

        double operator()(const std::vector<double>& path) const override {
            double logSum = 0.0;
            for (size_t i = 1; i < path.size(); ++i)
                logSum += std::log(path[i]);

            double geomMean = std::exp(logSum / (path.size() - 1));
            return std::max(K_ - geomMean, 0.0);
        }

        double operator()(double ST) const override {
            return std::max(K_ - ST, 0.0); 
        }

    private:
        double K_;
};
