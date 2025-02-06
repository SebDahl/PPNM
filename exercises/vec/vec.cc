#include <iostream>
#include <cmath>
#include "vec.h"

// class for 3D euclidean vectors
    vec::vec(double x, double y, double z) : x(x), y(y), z(z) {}

    vec vec::operator+(const vec &v) {
        return vec(x + v.x, y + v.y, z + v.z);
    }

    vec operator-(const vec &v) const {
        return vec(x - v.x, y - v.y, z - v.z);
    }

    vec operator*(double s) const {
        return vec(x * s, y * s, z * s);
    }

    vec operator/(double s) const {
        return vec(x / s, y / s, z / s);
    }

    vec operator-() const {
        return vec(-x, -y, -z);
    }

    double dot(const vec &v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    vec cross(const vec &v) const {
        return vec(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    double norm() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    double norm2() const {
        return x * x + y * y + z * z;
    }

    vec normalized() const {
        return *this / norm();
    }

    vec &operator+=(const vec &v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    vec &operator-=(const vec &v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }

    vec &operator*=(double s) {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    vec &operator/=(double s) {
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }

    bool operator==(const vec &v) const {
        return x == v.x && y == v.y && z == v.z;
}


std::ostream &operator<<(std::ostream &os, const vec &v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}
