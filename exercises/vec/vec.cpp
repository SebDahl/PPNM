#include <iostream>
#include <cmath>
#include <string>
#include "vec.h"

// Constructor: initializes the vector components.
vec::vec(double x, double y, double z) : x(x), y(y), z(z) {}

// In-place addition.
vec& vec::operator+=(const vec& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

// In-place subtraction.
vec& vec::operator-=(const vec& other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

// In-place scalar multiplication.
vec& vec::operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

// In-place scalar division.
vec& vec::operator/=(double s) {
    // Optionally: check for division by zero.
    x /= s;
    y /= s;
    z /= s;
    return *this;
}

// Dot product (member function).
double vec::dot(const vec& other) const {
    return x * other.x + y * other.y + z * other.z;
}

// Cross product. (also member function)
vec vec::cross(const vec& other) const {
    return vec( y * other.z - z * other.y,
                z * other.x - x * other.z,
                x * other.y - y * other.x );
}

// Print function for debugging.
void vec::print(std::string s) const {
    if (!s.empty())
        std::cout << s << ": ";
    std::cout << "{ " << x << ", " << y << ", " << z << " }" << std::endl;
}

// Friend: output stream operator.
std::ostream& operator<<(std::ostream& os, const vec& v) {
    os << "{ " << v.x << ", " << v.y << ", " << v.z << " }";
    return os;
}

// Non-member operator: vector addition.
vec operator+(const vec& v1, const vec& v2) {
    return vec(v1.x, v1.y, v1.z) += v2;
}

// Non-member operator: vector subtraction.
vec operator-(const vec& v1, const vec& v2) {
    return vec(v1.x, v1.y, v1.z) -= v2;
}

// Non-member operator: unary minus.
vec operator-(const vec& v) {
    return vec(-v.x, -v.y, -v.z);
}

// Non-member operator: vector multiplied by scalar.
vec operator*(const vec& v, double s) {
    return vec(v.x, v.y, v.z) *= s;
}

// Non-member operator: scalar multiplied by vector.
vec operator*(double s, const vec& v) {
    return vec(v.x, v.y, v.z) *= s;
}

// Non-member operator: vector divided by scalar.
vec operator/(const vec& v, double s) {
    return vec(v.x, v.y, v.z) /= s;
}

// Approximate equality check.
bool approx(const vec& v1, const vec& v2, double acc, double /*eps*/) {
    return (std::abs(v1.x - v2.x) <= acc) &&
           (std::abs(v1.y - v2.y) <= acc) &&
           (std::abs(v1.z - v2.z) <= acc);
}
