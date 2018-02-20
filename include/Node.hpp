#ifndef Node_hpp
#define Node_hpp

#include "Data.hpp"
#include "Vector.hpp"

class Mesh;

class Node {
public:
    long ID;

    Data data;
    Vector normal;

    Node(Mesh &mesh, double x, double y);
};

#endif