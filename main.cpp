#include <iostream>
#include <vector>
#include "particle.h"

int main() {
    using FissionSplit = Particle::FissionSplit;

    std::cout << "Testing fission_splits for masses 2 to 1000 with energy input = energy(mass=1)\n";

    // Compute energy for mass = 1
    Particle tiny(V2D(0,0));
    tiny.set_mass(1.f);
    float energy_input = tiny.energy();

    for (int mass = 2; mass <= 1000; mass++) {
        Particle p(V2D(0,0));
        p.set_mass(static_cast<float>(mass));

        try {
            auto splits = p.fission_splits(energy_input);

            std::cout << "Mass = " << mass << ", splits = " << splits.size() << "\n";

            for (const auto& split : splits) {
                auto [A1,A2,nu,prob] = split;
                std::cout << "  Split: A1=" << A1
                          << " A2=" << A2
                          << " neutrons=" << nu
                          << " probability=" << prob << "\n";
            }
        } catch (const std::exception &e) {
            std::cout << "Mass = " << mass << " threw exception: " << e.what() << "\n";
        }
    }

    std::cout << "Fission_splits test completed.\n";
    return 0;
}
