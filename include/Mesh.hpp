#ifndef Mesh_hpp
#define Mesh_hpp

#include "Node.hpp"
#include "Edge.hpp"
#include "Cell.hpp"

#include <vector>

class Mesh {
public:
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    std::vector<Cell> cells;
    std::map<std::pair<int, int>, int> mapNodesWithEdge;

    inline long getNewNodeID(){ return nodes.size();};
    int InitMesh(Mesh *mesh);
};

#endif
