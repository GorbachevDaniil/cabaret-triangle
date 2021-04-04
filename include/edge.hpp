#ifndef Edge_hpp
#define Edge_hpp

#include <vector>

#include "vector.hpp"

class Mesh;

class Edge {
   public:
    long ID;
    double length;
    bool boundEdge;

    Vector normal;

    std::vector<long> endNodeIDs;
    std::vector<long> nodeIDs;
    std::vector<long> innerNodeIDs;
    std::vector<long> usedNodeIDs;
    std::vector<long> cellIDs;

    Edge(){};
    Edge(Mesh &mesh, long ID, long startNodeID, long endNodeID, bool boundEdge, int innerNodeNum);

    long getAnotherEndNode(long ID);
    long getNearInnerNode(long ID);
    long getFarInnerNode(long ID);
};

#endif