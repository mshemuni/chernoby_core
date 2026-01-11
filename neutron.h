#pragma once
#ifndef CHERNOBYL_CORE_ATOM_H
#define CHERNOBYL_CORE_ATOM_H

#include "particle.h"

class Neutron : public Particle {
public:

    int mass = 1;
    float decay_probability = 0.f;

    Particle& set_mass(int v) noexcept {
        mass = 1;
        return *this;
    }

    Neutron& set_decay_probability(float v) noexcept {
        decay_probability = 0.f;
        return *this;
    }

    float energy() const noexcept {
        return ENERGY_MULTIPLIER * mass * velocity.dot(velocity);
    }

    std::vector<FissionSplit> bohr_wheeler(std::optional<float> energy_input = std::nullopt) const {
        return {
            FissionSplit{ mass, 0, 0, 1.0f }
        };
    }

};

#endif // CHERNOBYL_CORE_ATOM_H
