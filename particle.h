#pragma once
#ifndef CHERNOBYL_CORE_PARTICLE_H
#define CHERNOBYL_CORE_PARTICLE_H

#include "utils.h"
#include "v2d.h"
#include <cmath>
#include <random>
#include <optional>
#include <vector>
#include <tuple>
#include <algorithm>

class Particle {
public:
    V2D position;
    V2D velocity;
    V2D acceleration;

    static constexpr float ENERGY_POWER = 0.7f;
    static constexpr float ENERGY_MULTIPLIER = 10.f;

    float time_lived = 0.f;
    float time_to_live = 1.f;

    int mass = 1;
    float sigma = 1.f;

    float attraction_strength = 0.f;
    float absorption_probability = 0.1f;
    float decay_probability = 0.f;

    float received_energy = 0.f;
    bool unstable = false;

    float attraction_softening = 0.01f;
    float restitution = 0.f;
    float resistance = 0.f;
    float maximum_velocity = std::numeric_limits<float>::max();

    float fusion_probability = 0.f;

    // --- Construction ---
    explicit Particle(
        const V2D& position,
        V2D velocity = V2D::random(),
        V2D acceleration = V2D()
    ) noexcept
        : position(position),
          velocity(velocity),
          acceleration(acceleration)
    {}

    // ============================================================
    // Setters (fluent, validated)
    // ============================================================

    Particle& set_position(const V2D& p) noexcept {
        position = p;
        return *this;
    }

    Particle& set_velocity(const V2D& v) noexcept {
        velocity = v;
        return *this;
    }

    Particle& set_acceleration(const V2D& a) noexcept {
        acceleration = a;
        return *this;
    }

    Particle& set_time_lived(float v) noexcept {
        time_lived = (v < 0.f) ? 0.f : v;
        return *this;
    }

    Particle& set_time_to_live(float v) noexcept {
        if (v > 0.f) time_to_live = v;
        return *this;
    }

    Particle& set_mass(int v) noexcept {
        if (v > 0) mass = v;
        return *this;
    }


    Particle& set_sigma(float v) noexcept {
        if (v > 0.f) sigma = v;
        return *this;
    }

    Particle& set_attraction_strength(float v) noexcept {
        if (v >= 0.f) attraction_strength = v;
        return *this;
    }

    Particle& set_absorption_probability(float v) noexcept {
        absorption_probability = Utils::clamp(v, 0.0f, 1.0f);
        return *this;
    }

    Particle& set_decay_probability(float v) noexcept {
        decay_probability = Utils::clamp(v, 0.0f, 1.0f);
        return *this;
    }

    Particle& set_unstable(bool v) noexcept {
        unstable = v;
        return *this;
    }

    Particle& set_attraction_softening(bool v) noexcept {
        attraction_softening = v;
        return *this;
    }

    Particle& set_restitution(bool v) noexcept {
        restitution = v;
        return *this;
    }

    Particle& set_resistance(float v) noexcept {
        resistance =  Utils::clamp(v, 0.0f, 1.0f);
        return *this;
    }

    Particle& set_maximum_velocity(float v) noexcept {
        maximum_velocity = Utils::clamp(v, 0.0f, 1.0f);
        return *this;
    }

    Particle& set_fusion_probability(float v) noexcept {
        fusion_probability = Utils::clamp(v, 0.0f, 1.0f);
        return *this;
    }

    // ============================================================
    // Derived quantities
    // ============================================================

    float radius() const noexcept {
        return Utils::mass_to_radius(mass, sigma);
    }

    float kinetic_energy() const noexcept {
        return ENERGY_MULTIPLIER * Utils::kinetic_energy(mass, velocity);
    }

    float mass_energy() const noexcept {
        return ENERGY_MULTIPLIER * Utils::mass_energy(mass, ENERGY_POWER);
    }

    float stiffness() const noexcept {
        return 1.f / mass_energy();
    }

    // ============================================================
    // Lifecycle
    // ============================================================

    void age(float dt = 1.0f/60.0f) noexcept {
        time_lived += dt;
    }

    void receive_energy(float v) noexcept {
        received_energy += v;
    }

    bool is_dead() const noexcept {
        return time_lived >= time_to_live;
    }

    bool is_outside(int width, int height) const noexcept {
        return position.x <= 0.f || position.x >= width ||
               position.y <= 0.f || position.y >= height;
    }

    // ============================================================
    // Motion / integration
    // ============================================================

    void integrate(float dt = 1.0f/60.0f) noexcept {
        apply_drag();
        velocity = velocity.add(acceleration.scale(dt));
        clamp_velocity();
        position = position.add(velocity.scale(dt));
        reset_forces();
    }

    void clamp_velocity() noexcept {
        float m = velocity.magnitude();
        if (m > maximum_velocity)
            velocity = velocity.unit().scale(maximum_velocity);
    }

    void reset_forces() noexcept {
        acceleration = V2D();
    }

    // ============================================================
    // Interaction
    // ============================================================

    bool is_collided(const Particle& other) const noexcept {
        float d = position.distance(other.position);
        return d <= radius() + other.radius();
    }

    void apply_force(const V2D& force) noexcept {
        acceleration = acceleration.add(force.scale(1.f / mass));
    }

    void apply_drag() noexcept {
        velocity = velocity.scale(1.f - resistance);
    }

    void attract_to(const Particle& other) noexcept {
        const V2D direction = other.position.subtract(position);
        const float distance = direction.magnitude();

        if (distance <= 0.f)
            return;

        const float force_mag =
            attraction_strength *
            static_cast<float>(mass) *
            static_cast<float>(other.mass) /
            (distance * distance + attraction_softening);

        apply_force(direction.unit().scale(force_mag));
    }

    void bounce(Particle& other) noexcept {
        if (!is_collided(other))
            return;

        if (mass <= 0 || other.mass <= 0)
            return;

        const V2D normal = position.subtract(other.position).unit();
        const V2D rel_vel = velocity.subtract(other.velocity);

        const float vel_n = rel_vel.dot(normal);
        if (vel_n > 0.f)
            return;

        const float inv_m1 = 1.f / static_cast<float>(mass);
        const float inv_m2 = 1.f / static_cast<float>(other.mass);

        float j =
            -(1.f + (restitution + other.restitution) * 0.5f) * vel_n /
            (inv_m1 + inv_m2);

        const V2D impulse = normal.scale(j);

        velocity = velocity.add(impulse.scale(inv_m1));
        other.velocity = other.velocity.subtract(impulse.scale(inv_m2));
    }

    bool should_absorb(const float dt = 1.0f/60.0f) const noexcept {
        return Utils::time_scaled_bernoulli(dt, absorption_probability);
    }

    bool should_decay(const float dt = 1.0f/60.0f) const noexcept {
        return Utils::time_scaled_bernoulli(dt, decay_probability);
    }

    bool can_fuse(const Particle& other) const noexcept
    {
        if (!is_collided(other))
            return false;

        const float m1 = static_cast<float>(mass);
        const float m2 = static_cast<float>(other.mass);

        if (m1 <= 0.f || m2 <= 0.f)
            return false;

        const V2D rel_vel = velocity.subtract(other.velocity);

        V2D normal = position.subtract(other.position);
        float dist = normal.magnitude();

        if (dist < 1e-6f)
            return false;

        normal = normal.unit();

        const float v_n = rel_vel.dot(normal);

        if (v_n >= 0.f)
            return false;

        const float reduced_mass = (m1 * m2) / (m1 + m2);

        const float kinetic_energy =  random().set_mass(reduced_mass).set_velocity(rel_vel).kinetic_energy();

        const float available_energy = mass_energy() + other.mass_energy() + kinetic_energy;

        const float fused_mass_energy = random().set_mass(mass + other.mass). mass_energy();

        return available_energy >= fused_mass_energy;
    }

    bool should_fuse(const Particle& other, float dt = 1.0f / 60.0f) const noexcept
    {
        if (!can_fuse(other))
            return false;

        return Utils::time_scaled_bernoulli(
            dt,
            (fusion_probability + other.absorption_probability) * 0.5f
            );
    }


    Particle fuse_with(const Particle& other) const
    {
        if (!can_fuse(other)) {
            throw std::logic_error("Particle::fuse_with called when fusion is not possible");
        }

        const float m1 = static_cast<float>(mass);
        const float m2 = static_cast<float>(other.mass);
        const float total_mass = m1 + m2;

        const V2D fused_position =
            position.scale(m1)
            .add(other.position.scale(m2))
            .scale(1.f / total_mass);

        const V2D fused_velocity =
            velocity.scale(m1)
            .add(other.velocity.scale(m2))
            .scale(1.f / total_mass);

        Particle fused(fused_position);
        fused
            .set_velocity(fused_velocity)
            .set_mass(static_cast<int>(total_mass))
            .set_sigma((sigma + other.sigma) * 0.5f)
            .set_restitution((restitution + other.restitution) * 0.5f)
            .set_resistance((resistance + other.resistance) * 0.5f)
            .set_attraction_strength((attraction_strength + other.attraction_strength) * 0.5f)
            .set_fusion_probability((fusion_probability + other.fusion_probability) * 0.5f)
            .set_unstable(unstable && other.unstable);


        return fused;
    }

    // ============================================================
    // Fission (mass-based)
    // ============================================================

    using FissionSplit = std::tuple<int,int,int,float>;

    std::vector<FissionSplit> bohr_wheeler(std::optional<float> energy_input = std::nullopt) const
    {
        if (mass < 3) {
            return {
                FissionSplit{ mass, 0, 0, 1.0f }
            };
        }

        float deltaA0 = 0.18f * std::pow(mass, 2.f / 3.f);
        if (energy_input.has_value()) {
            const float symmetry_factor =
                std::max(0.6f, 1.f - energy_input.value() / 100.f);
            deltaA0 *= symmetry_factor;
        }

        float base_nu = 2.5f;
        if (energy_input.has_value())
            base_nu += 0.02f * energy_input.value();

        const float max_delta = deltaA0 * 2.f;
        const float step = std::max(1.f, deltaA0 / 3.f);

        std::vector<FissionSplit> results;
        std::vector<float> weights;

        for (float deltaA = -max_delta; deltaA <= max_delta; deltaA += step) {
            int nu = static_cast<int>(
                std::round(base_nu + std::abs(deltaA) * 0.1f)
            );

            int A1 = static_cast<int>(
                std::round((mass - nu) / 2.f + deltaA)
            );
            int A2 = (mass - nu) - A1;

            if (A1 <= 0 || A2 <= 0)
                continue;

            float weight =
                std::exp(-(deltaA * deltaA) /
                         (2.f * deltaA0 * deltaA0));

            results.emplace_back(A1, A2, nu, weight);
            weights.push_back(weight);
        }

        float total_weight = 0.f;
        for (float w : weights) total_weight += w;

        std::vector<FissionSplit> normalized;
        for (size_t i = 0; i < results.size(); ++i) {
            auto [A1, A2, nu, _] = results[i];
            normalized.emplace_back(
                A1, A2, nu, weights[i] / total_weight
            );
        }

        return normalized;
    }

    FissionSplit bohr_wheeler_split(std::optional<float> energy_input = std::nullopt) const {
        auto splits = bohr_wheeler(energy_input);

        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dist(0.f, 1.f);

        float r = dist(gen);
        float cum_prob = 0.f;

        for (const auto& s : splits) {
            cum_prob += std::get<3>(s);
            if (r <= cum_prob)
                return s;
        }

        return splits.back();
    }

    // ============================================================
    // Factory helpers
    // ============================================================

    static Particle random()  noexcept{
        return Particle(V2D::random());
    }

    static Particle from_values(float x, float y)  noexcept{
        return Particle(V2D(x, y));
    }
};

#endif // CHERNOBYL_CORE_PARTICLE_H
