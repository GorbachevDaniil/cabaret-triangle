#ifndef Node_hpp
#define Node_hpp

#include "Data.hpp"

class Mesh;

class Node {
   public:
    long ID;
    bool used;

    Data data;

    Node(){};
    Node(Mesh &mesh, double x, double y, bool used);
};

#endif
