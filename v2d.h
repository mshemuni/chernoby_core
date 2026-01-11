//
// Created by niaei on 9.01.2026.
//

#pragma once
#ifndef CHERNOBY_CORE_V2D_H
#define CHERNOBY_CORE_V2D_H

#include <algorithm>
#include <cmath>
#include <tuple>
#include <stdexcept>
#include <random>

class V2D {
public:
    static constexpr float PI = 3.14159265358979323846f;

    float x, y;

    constexpr V2D(float x = 0.f, float y = 0.f) noexcept
    : x(x), y(y) {}

    constexpr V2D add(const V2D& other) const noexcept {
        return V2D(x + other.x, y + other.y);
    }

    constexpr V2D subtract(const V2D& other) const noexcept {
        return V2D(x - other.x, y - other.y);
    }

    constexpr float dot(const V2D& other) const noexcept {
        return x * other.x + y * other.y;
    }

    constexpr V2D scale(float scalar) const noexcept {
        return V2D(x * scalar, y * scalar);
    }

    constexpr V2D multiply(float scalar) const noexcept {
        return scale(scalar);
    }
    
    V2D divide(float scalar) const {
        if (scalar == 0.f)
            throw std::invalid_argument("Cannot divide by zero");
        return scale(1.f / scalar);
    }

    float magnitude() const noexcept {
        return std::hypot(x, y);
    }

    float distance(const V2D& other = V2D()) const noexcept {
        return std::hypot(x - other.x, y - other.y);
    }

    V2D unit() const {
        float m = magnitude();
        if (m == 0.f)
            throw std::runtime_error("Cannot get unit of zero vector");
        return scale(1.f / m);
    }

    float angle_between(const V2D& other) const {
        float m1 = magnitude();
        float m2 = other.magnitude();
        if (m1 == 0.f || m2 == 0.f)
            throw std::runtime_error("Cannot compute angle with zero vector");

        float cos_angle = dot(other) / (m1 * m2);
        cos_angle = std::clamp(cos_angle, -1.f, 1.f);
        return std::acos(cos_angle) * 180.f / PI;
    }

    bool is_parallel(const V2D& other, float tolerance = 1e-6f) const noexcept {
        return std::fabs(x * other.y - y * other.x) < tolerance;
    }

    bool is_perpendicular(const V2D& other, float tolerance = 1e-6f) const noexcept {
        return std::fabs(dot(other)) < tolerance;
    }

    bool is_non_parallel(const V2D& other) const noexcept {
        return !(is_parallel(other) || is_perpendicular(other));
    }

    bool is_same(const V2D& other, float tolerance = 1e-6f) const noexcept {
        return distance(other) < tolerance;
    }

    V2D rotate(float angle_deg) const noexcept {
        float rad = angle_deg * PI / 180.f;
        float cos_a = std::cos(rad);
        float sin_a = std::sin(rad);
        return V2D(
            x * cos_a - y * sin_a,
            x * sin_a + y * cos_a
        );
    }

    V2D rotate(float angle_deg, const V2D& other) const noexcept {
        return subtract(other)
            .rotate(angle_deg)
            .add(other);
    }

    std::tuple<float, float> as_tuple(bool whole = false) const noexcept {
        if (whole)
            return { static_cast<int>(x), static_cast<int>(y) };
        return { x, y };
    }

    static V2D from_polar(float magnitude, float angle_deg) noexcept {
        float rad = angle_deg * PI / 180.f;
        return V2D(
            magnitude * std::cos(rad),
            magnitude * std::sin(rad)
        );
    }

    static V2D random() noexcept {
        static thread_local std::mt19937 rng{ std::random_device{}() };
        static std::uniform_real_distribution<float> angle_dist(0.f, 360.f);

        float angle = angle_dist(rng);
        return from_polar(1.f, angle);
    }


};

#endif //CHERNOBY_CORE_V2D_H