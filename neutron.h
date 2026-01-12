#pragma once

#include "particle.h"

class Neutron : public Particle {
public:

    // Mass and Radius
    const int mass = 1;

    // Nuclear Operations
    const float absorption_probability = 0.0f;
    const float decay_probability = 0.f;
    const float fusion_probability = 0.f;
    const bool unstable = false;

    // Classical Physics
    const float force_softening = 0.0f;
    const float restitution = 0.f;
    const float resistance = 0.f;
    float maximum_velocity = std::numeric_limits<float>::max();

    // --- Construction ---
    explicit Neutron(
        const V2D& position,
        const V2D velocity = V2D::random_vector(),
        const V2D acceleration = V2D()
        ) noexcept
        : Particle(position, velocity, acceleration) {}
    // ============================================================
    // Setters (fluent, validated), overwrite
    // ============================================================

    Neutron& set_mass(int value) noexcept override {
        return *this;
    }

    Neutron& set_absorption_probability(float value) noexcept override {
        return *this;
    }

    Neutron& set_decay_probability(float value) noexcept override {
        return *this;
    }

    Neutron& set_fusion_probability(float value) noexcept override {
        return *this;
    }

    Neutron& set_unstable(const bool value) noexcept override {
        return *this;
    }

    Neutron& set_force_softening(float value) noexcept override {
        return *this;
    }

    Neutron& set_restitution(float value) noexcept override {
        return *this;
    }

    Neutron& set_resistance(float value) noexcept override {
        return *this;
    }

    Neutron& set_maximum_velocity(float value) noexcept override {
        return *this;
    }

    // ============================================================
    // Derived quantities, overwrite
    // ============================================================

    [[nodiscard]]
    float radius() const noexcept override {
        return Utils::mass_to_radius(mass, sigma) * 0.5f;
    }

    [[nodiscard]]
    float mass_energy() const noexcept override {
        return 0.f;
    }

    // ============================================================
    // Lifecycle, overwrite
    // ============================================================

    void receive_energy(const float energy) noexcept override {}

    // ============================================================
    // Motion / integration, overwrite
    // ============================================================

    void clamp_velocity() noexcept override {}
    void apply_drag() noexcept override {}

    // ============================================================
    // Interaction, overwrite
    // ============================================================

    std::vector<FissionSplit> bohr_wheeler(const float energy_input) const override {
        return {
            FissionSplit{ mass, 0, 0, 1.0f }
        };
    }

    bool should_absorb(float dt = 1.0f/60.0f) const noexcept override {
        return false;
    }

    bool should_decay(float dt = 1.0f/60.0f) const noexcept override {
        return false;
    }

    [[nodiscard]]
    bool can_fuse(const Particle& other) const noexcept override {
        return false;
    }

    [[nodiscard]]
    bool should_fuse(const Particle& other, const float delta_time = 1.0f / 60.0f) const noexcept override {
        return false;
    }

};
