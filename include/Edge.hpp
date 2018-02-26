#ifndef Edge_hpp
#define Edge_hpp

#include <vector>

class Mesh;

class Edge {
public:
    long ID;
    long centerNodeID;
    double length;
    bool boundEdge;

    std::vector<long> nodeIDs;
    std::vector<long> cellIDs;

    Edge() {};
    Edge(Mesh &mesh, long id, long start_node_id, long end_node_id, bool boundEdge);
};

#endif