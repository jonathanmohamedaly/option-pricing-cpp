#include "../include/RandomGenerator.hpp"

RandomGenerator::RandomGenerator(unsigned int seed)
    : engine_(seed), dist_(0.0, 1.0) {}

double RandomGenerator::gaussian() {
    return dist_(engine_);
}