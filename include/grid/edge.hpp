#ifndef Edge_hpp
#define Edge_hpp

#include <vector>

#include "vector.hpp"

class Mesh;

class Edge {
public:
    long id;
    double length;
    bool is_bound;

    std::vector<long> end_node_ids;
    std::vector<long> node_ids;
    std::vector<long> inner_node_ids;
    std::vector<long> used_node_ids;
    std::vector<long> cell_ids;

    Edge(Mesh &mesh, long ID, long startNodeID, long endNodeID, bool boundEdge, int innerNodeNum);

    long getAnotherEndNode(long ID);
    long getNearInnerNode(long ID);
    long getFarInnerNode(long ID);
};

#endif