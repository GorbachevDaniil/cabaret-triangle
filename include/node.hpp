#ifndef Node_hpp
#define Node_hpp

#include <set>

#include "data.hpp"

class Mesh;

class Node {
public:
    long ID;
    bool used;
    bool boundNode;
    bool cellCenterNode;
    bool isApex;

    Data data;

    std::set<long> cellIDs;
    std::set<long> edgeIDs;

    Node() {};
    Node(Mesh &mesh, double x, double y, bool used, bool boundNode, bool cellCenterNode,
         bool isApex);
};

#endif
