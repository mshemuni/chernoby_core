#pragma once

#include <cmath>
#include <random>
#include <vector>
#include <tuple>
#include <algorithm>
#include <limits>

#include "utils.h"
#include "v2d.h"

class Particle {
public:

    V2D position;
    V2D velocity;
    V2D acceleration;

    // Energy Calculation
    static constexpr float ENERGY_POWER = 0.7f;
    static constexpr float ENERGY_MULTIPLIER = 10.f;

    // Timing
    float time_lived = 0.f;
    float time_to_live = 1.f;

    // Mass and Radius
    int mass = 1;
    float sigma = 1.f;

    // Nuclear Operations
    float absorption_probability = 0.1f;
    float decay_probability = 0.f;
    float fusion_probability = 0.f;
    float received_energy = 0.f;
    bool unstable = false;

    // Classical Physics
    float attraction_strength = 0.f;
    float force_softening = 0.01f;
    float restitution = 0.f;
    float resistance = 0.f;
    float maximum_velocity = std::numeric_limits<float>::max();

    // --- Construction ---
    virtual ~Particle() = default;

    explicit Particle(
        const V2D& position,
        const V2D velocity = V2D::random_vector(),
        const V2D acceleration = V2D()
    ) noexcept
        : position(position),
          velocity(velocity),
          acceleration(acceleration)
    {}

    // ============================================================
    // Setters (fluent, validated)
    // ============================================================

    virtual Particle& set_position(const V2D& _position) noexcept {
        position = _position;
        return *this;
    }

    virtual Particle& set_velocity(const V2D& _velocity) noexcept {
        velocity = _velocity;
        return *this;
    }

    virtual Particle& set_acceleration(const V2D& _acceleration) noexcept {
        acceleration = _acceleration;
        return *this;
    }

    virtual Particle& set_time_lived(const float value) noexcept {
        time_lived = (value < 0.f) ? 0.f : value;
        return *this;
    }

    virtual Particle& set_time_to_live(const float value) noexcept {
        if (value > 0.f) time_to_live = value;
        return *this;
    }

    virtual Particle& set_mass(const int value) noexcept {
        if (value > 0) mass = value;
        return *this;
    }

    virtual Particle& set_sigma(const float value) noexcept {
        if (value > 0.f) sigma = value;
        return *this;
    }

    virtual Particle& set_attraction_strength(const float value) noexcept {
        attraction_strength = value;
        return *this;
    }

    virtual Particle& set_absorption_probability(const float value) noexcept {
        absorption_probability = Utils::clamp(value, 0.0f, 1.0f);
        return *this;
    }

    virtual Particle& set_decay_probability(const float value) noexcept {
        decay_probability = Utils::clamp(value, 0.0f, 1.0f);
        return *this;
    }

    virtual Particle& set_unstable(const bool value) noexcept {
        unstable = value;
        return *this;
    }

    virtual Particle& set_force_softening(const float value) noexcept {
        force_softening = value;
        return *this;
    }

    virtual Particle& set_restitution(const float value) noexcept {
        restitution = value;
        return *this;
    }

    virtual Particle& set_resistance(const float value) noexcept {
        resistance =  Utils::clamp(value, 0.0f, 1.0f);
        return *this;
    }

    virtual Particle& set_maximum_velocity(const float value) noexcept {
        maximum_velocity = std::max(0.f, value);
        return *this;
    }

    virtual Particle& set_fusion_probability(const float value) noexcept {
        fusion_probability = Utils::clamp(value, 0.0f, 1.0f);
        return *this;
    }

    // ============================================================
    // Derived quantities
    // ============================================================

    [[nodiscard]]
    virtual float radius() const noexcept {
        return Utils::mass_to_radius(mass, sigma);
    }

    [[nodiscard]]
    virtual float kinetic_energy() const noexcept {
        return ENERGY_MULTIPLIER * Utils::kinetic_energy(mass, velocity);
    }

    [[nodiscard]]
    virtual float mass_energy() const noexcept {
        return ENERGY_MULTIPLIER * Utils::mass_energy(mass, ENERGY_POWER);
    }

    // ============================================================
    // Lifecycle
    // ============================================================

    virtual void age(const float delta_time) noexcept {
        time_lived += delta_time;
    }

    virtual void age() noexcept {
        age(1.0f/60.0f);
    }

    virtual void receive_energy(const float energy) noexcept {
        received_energy += energy;
    }

    [[nodiscard]]
    virtual bool is_dead() const noexcept {
        return time_lived >= time_to_live;
    }

    [[nodiscard]]
    virtual bool is_outside(const int width, const int height) const noexcept {
        return position.x <= 0.f || position.x >= static_cast<float>(width) ||
               position.y <= 0.f || position.y >= static_cast<float>(height);
    }

    // ============================================================
    // Motion / integration
    // ============================================================

    virtual void integrate(const float delta_time) noexcept {
        apply_drag();
        velocity = velocity.add(acceleration.scale(delta_time));
        clamp_velocity();
        position = position.add(velocity.scale(delta_time));
        reset_forces();
    }

    virtual void integrate() noexcept {
        integrate(1.0f/60.0f);
    }

    virtual void clamp_velocity() noexcept {
        float _mass = velocity.magnitude();
        if (_mass > maximum_velocity)
            velocity = velocity.unit().scale(maximum_velocity);
    }

    virtual void reset_forces() noexcept {
        acceleration = V2D();
    }

    // ============================================================
    // Interaction
    // ============================================================

    [[nodiscard]]
    virtual bool is_collided(const Particle& other) const noexcept {
        float d = position.distance(other.position);
        return d <= radius() + other.radius();
    }

    virtual void apply_force(const V2D& force) noexcept {
        acceleration = acceleration.add(
            force.scale(1.f / (static_cast<float>(mass) * force_softening))
        );
    }

    virtual void apply_drag() noexcept {
        velocity = velocity.scale(std::max(0.f, 1.f - resistance));
    }

    virtual void attract_to(const Particle& other) noexcept {
        const V2D direction = other.position.subtract(position);
        const float distance = direction.magnitude();

        if (distance <= 0.f)
            return;

        const float force_mag =
            attraction_strength *
            static_cast<float>(mass) *
            static_cast<float>(other.mass) /
            (distance * distance);

        apply_force(direction.unit().scale(force_mag));
    }

    virtual void bounce(Particle& other) noexcept {
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

    [[nodiscard]]
    virtual bool should_absorb(const float delta_time) const noexcept {
        return Utils::time_scaled_bernoulli(delta_time, absorption_probability);
    }

    [[nodiscard]]
    virtual bool should_absorb() const noexcept {
        return should_absorb(1.0f/60.0f);
    }

    [[nodiscard]]
    virtual bool should_decay(const float delta_time) const noexcept {
        return Utils::time_scaled_bernoulli(delta_time, decay_probability);
    }

    [[nodiscard]]
    virtual bool should_decay() const noexcept {
        return should_decay(1.0f/60.0f);
    }

    [[nodiscard]]
    virtual bool can_fuse(const Particle& other) const noexcept
    {
        if (!is_collided(other))
            return false;

        const auto m1 = static_cast<float>(mass);
        const auto m2 = static_cast<float>(other.mass);

        if (m1 <= 0.f || m2 <= 0.f)
            return false;

        const V2D rel_vel = velocity.subtract(other.velocity);

        V2D normal = position.subtract(other.position);
        const float dist = normal.magnitude();

        if (dist < 1e-6f)
            return false;

        normal = normal.unit();

        const float v_n = rel_vel.dot(normal);

        if (v_n >= 0.f)
            return false;

        const float reduced_mass = (m1 * m2) / (m1 + m2);

        const float kinetic_energy =  random_particle()
                                        .set_mass(static_cast<int>(reduced_mass))
                                        .set_velocity(rel_vel)
                                        .kinetic_energy();

        const float available_energy = mass_energy() + other.mass_energy() + kinetic_energy;

        const float fused_mass_energy = random_particle().set_mass(mass + other.mass). mass_energy();

        return available_energy >= fused_mass_energy;
    }

    [[nodiscard]]
    virtual bool should_fuse(const Particle& other, const float delta_time) const noexcept
    {
        if (!can_fuse(other))
            return false;

        return Utils::time_scaled_bernoulli(
            delta_time,
            (fusion_probability + other.absorption_probability) * 0.5f
            );
    }

    [[nodiscard]]
    virtual bool should_fuse(const Particle& other) const noexcept
    {
        return should_fuse(other, 1.0f/60.0f);
    }

    [[nodiscard]]
    virtual Particle fuse_with(const Particle& other) const
    {
        if (!can_fuse(other)) {
            throw std::logic_error("Particle::fuse_with called when fusion is not possible");
        }

        const auto m1 = static_cast<float>(mass);
        const auto m2 = static_cast<float>(other.mass);
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

    [[nodiscard]]
    virtual std::vector<FissionSplit> bohr_wheeler(const float energy_input) const {
        if (mass < 3) {
            return {
                FissionSplit{ mass, 0, 0, 1.0f }
            };
        }

        float deltaA0 = 0.18f * std::pow(static_cast<float>(mass), 2.f / 3.f);

        const float symmetry_factor =
                std::max(0.6f, 1.f - energy_input / 100.f);
        deltaA0 *= symmetry_factor;

        float base_nu = 2.5f + 0.02f * energy_input;

        const float max_delta = deltaA0 * 2.f;
        const float step = std::max(1.f, deltaA0 / 3.f);

        std::vector<FissionSplit> results;
        std::vector<float> weights;

        const int steps = static_cast<int>(std::floor((2 * max_delta) / step));
        for (int i = 0; i <= steps; ++i) {
            const float deltaA = -max_delta + static_cast<float>(i) * step;

            int nu = static_cast<int>(
                std::round(base_nu + std::abs(deltaA) * 0.1f)
            );

            int A1 = static_cast<int>(
                std::round((static_cast<float>(mass) - static_cast<float>(nu)) / 2.f + deltaA)
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
        for (const float w : weights) total_weight += w;

        std::vector<FissionSplit> normalized;
        for (size_t i = 0; i < results.size(); ++i) {
            auto [A1, A2, nu, _] = results[i];
            normalized.emplace_back(
                A1, A2, nu, weights[i] / total_weight
            );
        }

        return normalized;
    }

    [[nodiscard]]
    virtual FissionSplit bohr_wheeler_split(const float energy_input) const {
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

    static Particle random_particle()  noexcept{
        return Particle(V2D::random_vector());
    }

    static Particle from_values(const float x, const float y)  noexcept{
        return Particle(V2D(x, y));
    }
};
