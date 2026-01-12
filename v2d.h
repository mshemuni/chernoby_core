#pragma once

#include <algorithm>
#include <cmath>
#include <tuple>
#include <stdexcept>
#include <random>

class V2D {
public:
    static constexpr float PI = 3.14159265358979323846f;

    float x, y;

    explicit constexpr V2D(const float x = 0.f, const float y = 0.f) noexcept
    : x(x), y(y) {}

    [[nodiscard]]
    constexpr V2D add(const V2D& other) const noexcept {
        return V2D(x + other.x, y + other.y);
    }

    [[nodiscard]]
    constexpr V2D subtract(const V2D& other) const noexcept {
        return V2D(x - other.x, y - other.y);
    }

    [[nodiscard]]
    constexpr float dot(const V2D& other) const noexcept {
        return x * other.x + y * other.y;
    }

    [[nodiscard]]
    constexpr V2D scale(const float scalar) const noexcept {
        return V2D(x * scalar, y * scalar);
    }

    [[nodiscard]]
    constexpr V2D multiply(const float scalar) const noexcept {
        return scale(scalar);
    }

    [[nodiscard]]
    V2D divide(const float scalar) const {
        if (scalar == 0.f)
            throw std::invalid_argument("Cannot divide by zero");
        return scale(1.f / scalar);
    }

    [[nodiscard]]
    float magnitude() const noexcept {
        return std::hypot(x, y);
    }

    [[nodiscard]]
    float distance(const V2D& other = V2D()) const noexcept {
        return std::hypot(x - other.x, y - other.y);
    }

    [[nodiscard]]
    V2D unit() const {
        const float m = magnitude();
        if (m == 0.f)
            throw std::runtime_error("Cannot get unit of zero vector");
        return scale(1.f / m);
    }

    [[nodiscard]]
    float angle_between(const V2D& other) const {
        const float m1 = magnitude();
        const float m2 = other.magnitude();
        if (m1 == 0.f || m2 == 0.f)
            throw std::runtime_error("Cannot compute angle with zero vector");

        float cos_angle = dot(other) / (m1 * m2);
        cos_angle = std::clamp(cos_angle, -1.f, 1.f);
        return std::acos(cos_angle) * 180.f / PI;
    }

    [[nodiscard]]
    bool is_parallel(const V2D& other, const float tolerance = 1e-6f) const noexcept {
        return std::fabs(x * other.y - y * other.x) < tolerance;
    }

    [[nodiscard]]
    bool is_perpendicular(const V2D& other, const float tolerance = 1e-6f) const noexcept {
        return std::fabs(dot(other)) < tolerance;
    }

    [[nodiscard]]
    bool is_non_parallel(const V2D& other) const noexcept {
        return !(is_parallel(other) || is_perpendicular(other));
    }

    [[nodiscard]]
    bool is_same(const V2D& other, const float tolerance = 1e-6f) const noexcept {
        return distance(other) < tolerance;
    }

    [[nodiscard]]
    V2D rotate(const float angle_deg) const noexcept {
        const float rad = angle_deg * PI / 180.f;
        const float cos_a = std::cos(rad);
        const float sin_a = std::sin(rad);
        return V2D(
            x * cos_a - y * sin_a,
            x * sin_a + y * cos_a
        );
    }

    [[nodiscard]]
    V2D rotate(const float angle_deg, const V2D& other) const noexcept {
        return subtract(other)
            .rotate(angle_deg)
            .add(other);
    }

    [[nodiscard]]
    V2D whole() const noexcept {
        return V2D(std::trunc(x), std::trunc(y));
    }

    [[nodiscard]]
    std::tuple<float, float> as_tuple(const bool whole = false) const noexcept {
        if (whole)
            return { static_cast<int>(x), static_cast<int>(y) };
        return { x, y };
    }

    static V2D from_polar(const float magnitude, const float angle_deg) noexcept {
        const float rad = angle_deg * PI / 180.f;
        return V2D(
            magnitude * std::cos(rad),
            magnitude * std::sin(rad)
        );
    }

    static V2D random_vector() noexcept {
        static thread_local std::mt19937 rng{ std::random_device{}() };
        static std::uniform_real_distribution<float> angle_dist(0.f, 360.f);

        const float angle = angle_dist(rng);
        return from_polar(1.f, angle);
    }

};
