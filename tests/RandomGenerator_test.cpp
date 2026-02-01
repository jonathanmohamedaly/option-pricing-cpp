#include <iostream>
#include "../include/RandomGenerator.hpp"

int main() {
    RandomGenerator rng; 
    for (int i = 0; i < 5; ++i) {
        std::cout << "Gaussien " << i+1 << " : " << rng.gaussian() << std::endl;
    }
}