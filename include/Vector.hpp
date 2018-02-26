#ifndef Vector_hpp
#define Vector_hpp

#include <math.h>

class Vector {
public:
    double x;
    double y;

    inline void clear() { x = 0.;  y = 0.; }
    inline void set(Vector &c) { x = c.x; y = c.y; }
    inline void set(double *c) { x = c[0]; y = c[1]; }
    inline void set(double c0, double c1) {x = c0; y = c1; }
    inline void setMult(Vector &c, double d) { x = d * c.x; y = d * c.y; }
    inline void setPlus(Vector &c1, Vector &c2) { x=c1.x + c2.x; y=c1.y + c2.y; }
    inline void setMinus(Vector &c1, Vector &c2) { x=c1.x - c2.x; y=c1.y - c2.y; }
    inline void mult(double d){ x *= d; y *= d; }
    inline void plus(Vector &c){ x += c.x; y += c.y; }
    inline void minus(Vector &c){ x -= c.x; y -= c.y; }
    inline void plusMult(Vector &c, double d){ x += d * c.x; y += d * c.y; }
    inline void minusMult(Vector &c, double d){ x -= d * c.x; y -= d * c.y; }

    inline double scalar(Vector &c){ return (x * c.x + y * c.y); }
    inline double length(){ return sqrt(x * x + y * y); }
    inline double length2(){ return (x * x + y * y); }

    static double Scalar(Vector &v1, Vector &v2) { return v1.x * v2.x + v1.y * v2.y; }
    static double Scalar(double *d, Vector &v) { return d[0] * v.x + d[1] * v.y; }
    static double Scalar(double *d1, double *d2) { return d1[0] * d2[0] + d1[1] * d2[1]; }
    static double Length(double x, double y) { return sqrt(x * x + y * y); }
    static double GetDistance(Vector &v1, Vector &v2) { return Vector::Length(v2.x - v1.x,v2.y - v1.y); }
};

#endif