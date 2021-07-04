#ifndef Vector_hpp
#define Vector_hpp

#include <cmath>

class Vector {
public:
    double x;
    double y;

    Vector() : x(0), y(0) {};
    Vector(double x, double y) : x(x), y(y) {};

    double length() const { return sqrt(x * x + y * y); }
    static double length(double x, double y) { return sqrt(x * x + y * y); }
};

#endif
