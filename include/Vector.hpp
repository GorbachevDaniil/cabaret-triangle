#ifndef Vector_hpp
#define Vector_hpp

class Vector {
public:
    double x;
    double y;

    inline void Clear() { x = 0.;  y = 0.; }
    inline void Set(Vector &c) { x = c.x; y = c.y; }
    inline void Set(double *c) { x = c[0]; y = c[1]; }
    inline void Set(double c0, double c1) {x = c0; y = c1; }
    inline void SetMult(Vector &c, double d) { x = d * c.x; y = d * c.y; }
    inline void SetPlus(Vector &c1, Vector &c2) { x=c1.x + c2.x; y=c1.y + c2.y; }
    inline void SetMinus(Vector &c1, Vector &c2) { x=c1.x - c2.x; y=c1.y - c2.y; }
    inline void Mult(double d){ x *= d; y *= d; }
    inline void Plus(Vector &c){ x += c.x; y += c.y; }
    inline void Minus(Vector &c){ x -= c.x; y -= c.y; }
    inline void PlusMult(Vector &c, double d){ x += d * c.x; y += d * c.y; }
    inline void MinusMult(Vector &c, double d){ x -= d * c.x; y -= d * c.y; }

    inline double Scalar(Vector &c){ return (x * c.x + y * c.y); }
    inline double Length(){ return sqrt(x * x + y * y); }
    inline double Length2(){ return (x * x + y * y); }

    static double ScalarProduct(Vector &v1, Vector &v2) { return v1.x * v2.x + v1.y * v2.y; }
    static double ScalarProduct(double *d, Vector &v) { return d[0] * v.x + d[1] * v.y; }
    static double ScalarProduct(double *d1, double *d2) { return d1[0] * d2[0] + d1[1] * d2[1]; }
    static double Length(double x, double y, double z) { return sqrt(x * x + y * y); }
    static double GetDistance(Vector &v1, Vector &v2) { return Vector::Length(v2.x - v1.x,v2.y - v1.y); }
};

#endif