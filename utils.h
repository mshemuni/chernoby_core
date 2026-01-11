#ifndef CHERNOBY_CORE_UTILS_H
#define CHERNOBY_CORE_UTILS_H

#include <cmath>
#include <random>
#include <algorithm>

#include "v2d.h"

class Utils {
public:
    static float clamp(const float x, float min, float max) noexcept {
        if (min > max) std::swap(min, max);
        return std::clamp(x, min, max);
    }

    static bool time_scaled_bernoulli(const float time_scale,
                                      const float probability) noexcept
    {
        const float p = 1.0f - std::pow(1.0f - probability, time_scale);

        static thread_local std::mt19937 rng{ std::random_device{}() };
        static thread_local std::uniform_real_distribution<float> dist(0.0f, 1.0f);

        return dist(rng) < p;
    }

    static float mass_to_radius(const int mass, const float sigma = 1.0f) noexcept {
        return std::sqrt(
            static_cast<float>(mass) /
            (static_cast<float>(M_PI) * sigma)
        );
    }

    static float kinetic_energy(const int mass, const V2D velocity) noexcept {
        return mass * velocity.dot(velocity);
    }

    static float mass_energy(const int mass, const float energy_power=0.7f) noexcept {
        return std::pow(static_cast<float>(mass), energy_power);
    }
};

#endif // CHERNOBY_CORE_UTILS_H
