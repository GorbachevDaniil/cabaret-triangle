#ifndef Data_hpp
#define Data_hpp

#include "Vector.hpp"

class Data {
public:
    Vector coords;

    Vector u0;
    Vector u1;
    Vector u2;
    double phi0;
    double phi1;
    double phi2;
};

#endif