#pragma once

#include <algorithm>

#include "v2d.h"
#include "particle.h"

class Rod {
public:
    V2D start;
    V2D end;

    // Physical
    static constexpr float LOWER_STEP = 0.01f;
    float thickness = 15;

    // Nuclear Operations
    float absorption_probability = 0.1f;

    // ------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------

    constexpr Rod() noexcept = default;

    constexpr Rod(const V2D& start, const V2D& end) noexcept
        : start(start), end(end) {}

    // ------------------------------------------------------------
    // Setters
    // ------------------------------------------------------------

    Rod& set_start(const V2D& _start) noexcept {
        start = _start;
        return *this;
    }

    Rod& set_end(const V2D& _end) noexcept {
        end = _end;
        return *this;
    }

    Rod& set_thickness(const float _thickness) noexcept {
        if (_thickness>0) thickness=_thickness;
        return *this;
    }

    // ------------------------------------------------------------
    // Collision
    // ------------------------------------------------------------

    bool is_collided(const Particle& _particle) const noexcept {
        const V2D difference = end.subtract(start);
        const float len2 = difference.dot(difference);

        if (len2 <= 0.f) {
            return start.distance(_particle.position) <= _particle.radius();
        }

        float t = (_particle.position.subtract(start)).dot(difference) / len2;
        t = std::clamp(t, 0.f, 1.f);

        const V2D closest = start.add(difference.scale(t));

        return closest.distance(_particle.position) <= _particle.radius();
    }

    // ============================================================
    // Interaction
    // ============================================================

    [[nodiscard]]
    bool should_absorb(const float delta_time) const noexcept {
        return Utils::time_scaled_bernoulli(delta_time, absorption_probability);
    }

    [[nodiscard]]
    bool should_absorb() const noexcept {
        return should_absorb(1.0f/60.0f);
    }
};
