#ifndef Vector_hpp
#define Vector_hpp

#include <cmath>

class Vector {
   public:
    double x;
    double y;

    Vector() : x(0), y(0){};
    Vector(double x, double y) : x(x), y(y){};
    Vector(const Vector &v) {
        x = v.x;
        y = v.y;
    }

    double operator*(Vector const &other) { return x * other.x + y * other.y; }

    Vector operator*(double coef) {
        Vector v(x, y);
        v.x *= coef;
        v.y *= coef;
        return v;
    }

    Vector operator/(double coef) {
        Vector v(x, y);
        v.x /= coef;
        v.y /= coef;
        return v;
    }

    Vector operator+(Vector const &other) { return Vector(x + other.x, y + other.y); }

    Vector operator-(Vector const &other) { return Vector(x - other.x, y - other.y); }

    double length() { return sqrt(x * x + y * y); }
    static double length(double x, double y) { return sqrt(x * x + y * y); }
};

#endif
