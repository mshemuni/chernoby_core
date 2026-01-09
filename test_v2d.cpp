#include <cassert>
#include <cmath>
#include <iostream>

#include "v2d.h"

static constexpr float EPS = 1e-6f;

static bool feq(float a, float b, float eps = EPS) {
    return std::fabs(a - b) < eps;
}

int main() {

    // --- construction ---
    {
        V2D v;
        assert(feq(v.x, 0.f));
        assert(feq(v.y, 0.f));

        V2D u(3.f, 4.f);
        assert(feq(u.x, 3.f));
        assert(feq(u.y, 4.f));
    }

    // --- add / subtract ---
    {
        V2D a(5.f, 4.f);
        V2D b(2.f, 1.f);

        V2D c = a.add(b);
        assert(feq(c.x, 7.f));
        assert(feq(c.y, 5.f));

        V2D d = a.subtract(b);
        assert(feq(d.x, 3.f));
        assert(feq(d.y, 3.f));
    }

    // --- scale / multiply ---
    {
        V2D v(2.f, -3.f);
        V2D r1 = v.scale(2.f);
        V2D r2 = v.multiply(2.f);

        assert(feq(r1.x, 4.f));
        assert(feq(r1.y, -6.f));
        assert(feq(r2.x, 4.f));
        assert(feq(r2.y, -6.f));
    }

    // --- dot product ---
    {
        V2D a(1.f, 0.f);
        V2D b(0.f, 1.f);
        assert(feq(a.dot(b), 0.f));

        V2D c(1.f, 2.f);
        V2D d(3.f, 4.f);
        assert(feq(c.dot(d), 11.f)); // 1*3 + 2*4
    }

    // --- magnitude / distance ---
    {
        V2D v(3.f, 4.f);
        assert(feq(v.magnitude(), 5.f));

        V2D a(1.f, 1.f);
        V2D b(4.f, 5.f);
        assert(feq(a.distance(b), 5.f)); // sqrt(3^2 + 4^2)
    }

    // --- unit vector ---
    {
        V2D v(3.f, 4.f);
        V2D u = v.unit();
        assert(feq(u.magnitude(), 1.f));
    }

    // --- angle between (in degrees, convert if needed) ---
    {
        V2D a(1.f, 0.f);
        V2D b(0.f, 1.f);
        float angle = a.angle_between(b); // ensure this returns degrees
        assert(feq(angle, 90.f));
    }

    // --- parallel / perpendicular ---
    {
        V2D a(2.f, 2.f);
        V2D b(4.f, 4.f);
        assert(a.is_parallel(b));

        V2D c(1.f, 0.f);
        V2D d(0.f, 1.f);
        assert(c.is_perpendicular(d));
        assert(c.is_non_parallel(d) == false);
    }

    // --- rotate about origin (degrees) ---
    {
        V2D v(1.f, 0.f);
        V2D r = v.rotate(90.f); // assuming rotate() uses degrees
        assert(feq(r.x, 0.f));
        assert(feq(r.y, 1.f));
    }

    // --- rotate about arbitrary point ---
    {
        V2D p(3.f, 2.f);
        V2D center(1.f, 1.f);
        V2D r = p.rotate(90.f, center); // assuming degrees
        assert(feq(r.x, 0.f));
        assert(feq(r.y, 3.f));
    }

    // --- random unit vector ---
    {
        V2D v1 = V2D::random();
        V2D v2 = V2D::random();

        // magnitude should be 1
        assert(feq(v1.magnitude(), 1.f));
        assert(feq(v2.magnitude(), 1.f));

        // not zero vector
        assert(!(feq(v1.x, 0.f) && feq(v1.y, 0.f)));
        assert(!(feq(v2.x, 0.f) && feq(v2.y, 0.f)));

        // extremely unlikely to be identical
        if (v1.is_same(v2)) {
            std::cout << "Warning: random vectors identical (very rare)\n";
        }
    }

    std::cout << "All V2D tests passed\n";
    return 0;
}
