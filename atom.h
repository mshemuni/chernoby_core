#pragma once
#ifndef CHERNOBYL_CORE_ATOM_H
#define CHERNOBYL_CORE_ATOM_H

#include "particle.h"

class Moderator : public Particle {
public:
    float absorption_probability = 0.f;

    Moderator& set_absorption_probability(float v) noexcept {
        absorption_probability = 0;
        return *this;
    }

    bool is_dead() const noexcept {
        return false;
    }
};

#endif // CHERNOBYL_CORE_ATOM_H
