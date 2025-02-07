#include <iostream>
#include <cmath>
#include "vec.h"

vec::vec(double x, double y, double z) : x(x), y(y), z(z) {}

vec vec::operator+ (const vec& v2) {
    return vec(x + v2.x, y + v2.y, z + v2.z);
}
vec vec::operator- (const vec& v2) {
    return vec(x - v2.x, y - v2.y, z - v2.z);
}
vec vec::operator* (double s) {
    return vec(x * s, y * s, z * s);
}

double vec::dot(const vec& v2) const {
    return x * v2.x + y * v2.y + z * v2.z;
}

vec vec::cross (const vec& v2) const {
    return vec(y * v2.z - z * v2.y,
               z * v2.x - x * v2.z,
               x * v2.y - y * v2.x);
}

vec vec::operator- () {
    return vec(-x, -y, -z);
}

vec& vec::operator=(const vec& v2) {
    if (this != &v2) {  // Self-assignment check
        x = v2.x;
        y = v2.y;
        z = v2.z;
    }
    return *this;
}
