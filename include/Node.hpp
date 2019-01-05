#ifndef Node_hpp
#define Node_hpp

#include "Data.hpp"

#include <set>

class Mesh;

class Node {
   public:
    long ID;
    bool used;
    bool boundNode;
    bool cellCenterNode;
    bool phase2Calculated;

    Data data;

    std::set<long> cellIDs;

    Node(){};
    Node(Mesh &mesh, double x, double y, bool used, bool boundNode, bool cellCenterNode);
};

#endif
