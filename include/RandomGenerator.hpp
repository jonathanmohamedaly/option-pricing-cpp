#pragma once
#include <random>


class RandomGenerator {

    public:

        RandomGenerator(unsigned int seed = 42);
        double gaussian();

    private:
    std::mt19937 engine_;                  
    std::normal_distribution<double> dist_; 
};