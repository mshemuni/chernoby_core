#pragma once

#include "particle.h"

class Atom : public Particle {
public:

    // Nuclear Operations
    float absorption_probability = 0.f;

    // Classical Physics
    float attraction_strength = 0.f;

    // --- Construction ---
    explicit Atom(
        const V2D& position,
        const V2D velocity = V2D::random_vector(),
        const V2D acceleration = V2D()
        ) noexcept
        : Particle(position, velocity, acceleration) {}

    // ============================================================
    // Setters (fluent, validated), Overwrite
    // ============================================================

    Atom& set_absorption_probability(const float value) noexcept override {
        absorption_probability = 0;
        return *this;
    }

    Atom& set_attraction_strength(const float value) noexcept override {
        attraction_strength = 0;
        return *this;
    }

    // ============================================================
    // Lifecycle, Overwrite
    // ============================================================

    bool is_dead() const noexcept override {
        return false;
    }
};
