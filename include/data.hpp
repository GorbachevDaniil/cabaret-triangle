#ifndef Data_hpp
#define Data_hpp

#include <vector>

#include "vector.hpp"

class Mesh;

class Data {
public:
    long nodeID;
    Vector coords;

    Mesh *mesh;

    double getS0(int pos);
    double getS1(int pos);
    double getS2(int pos);

    Vector getV0(int pos);
    Vector getV1(int pos);
    Vector getV2(int pos);
};

#endif
