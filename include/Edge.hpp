#ifndef Edge_hpp
#define Edge_hpp

#include <vector>

#include "Vector.hpp"

class Mesh;

class Edge {
   public:
    long ID;
    double length;
    bool boundEdge;

    Vector normal;

    std::vector<long> endNodeIDs;
    std::vector<long> nodeIDs;
    std::vector<long> cellIDs;

    Edge(){};
    Edge(Mesh &mesh, long ID, long startNodeID, long endNodeID, bool boundEdge, int innerNodeNum);

    std::vector<long> getUsedNodes(Mesh &mesh);
    std::vector<long> getInnerNodes();
};

#endif