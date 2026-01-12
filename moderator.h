#pragma once

#include "particle.h"

class Moderator : public Particle {
public:

    // Nuclear Operations
    const float absorption_probability = 1.f;
    const float decay_probability = 0.f;

    // --- Construction ---
    explicit Moderator(
        const V2D& position,
        const V2D velocity = V2D::random_vector(),
        const V2D acceleration = V2D()
        ) noexcept
        : Particle(position, velocity, acceleration) {}

    // ============================================================
    // Setters (fluent, validated)
    // ============================================================

    Moderator& set_absorption_probability(const float value) noexcept override {
        return *this;
    }

    Moderator& set_decay_probability(const float value) noexcept override {
        return *this;
    }

    // ============================================================
    // Derived quantities
    // ============================================================

    [[nodiscard]]
    float radius() const noexcept override {
        return 1.5f * Particle::radius();
    }

    [[nodiscard]]
    float mass_energy() const noexcept override {
        return 2.f * Particle::mass_energy();
    }
};
