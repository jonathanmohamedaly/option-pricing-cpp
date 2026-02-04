#pragma once
#include <random>


class RandomGenerator {

    public:

        RandomGenerator();
        double gaussian();

    private:
    std::mt19937 engine_;                  
    std::normal_distribution<double> dist_; 
};