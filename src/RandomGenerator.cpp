#include "../include/RandomGenerator.hpp"

RandomGenerator::RandomGenerator()
    : engine_(std::random_device{}()), dist_(0.0, 1.0) {}

double RandomGenerator::gaussian() {
    return dist_(engine_);
}