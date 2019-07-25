#ifndef Mesh_hpp
#define Mesh_hpp

#include "Cell.hpp"
#include "Edge.hpp"
#include "Node.hpp"

#include <vector>

class Mesh {
   public:
    int edgeInnerNodesNumber;
    bool apexNodesUsed;

    // list of scalar variables on n layer of time
    std::vector<std::vector<double>> s0;
    // list of scalar variables on n+1/2 layer of time
    std::vector<std::vector<double>> s1;
    // list of scalar variables on n+1 layer of time
    std::vector<std::vector<double>> s2;

    // list of vector variables on n layer of time
    std::vector<std::vector<Vector>> v0;
    // list of vector variables on n+1/2 layer of time
    std::vector<std::vector<Vector>> v1;
    // list of vector variables on n+1 layer of time
    std::vector<std::vector<Vector>> v2;

    std::vector<Node> nodes;
    std::vector<Edge> edges;
    std::vector<Cell> cells;

    std::map<std::pair<int, int>, int> mapNodesWithEdge;
    std::map<long, std::vector<double>> edgeIDToDivs;
    std::map<long, std::vector<Vector>> cellIDToGrads;

    Mesh(int edgeInnerNodesNumber, bool apexNodesUsed) {
        this->edgeInnerNodesNumber = edgeInnerNodesNumber;
        this->apexNodesUsed = apexNodesUsed;
    };

    inline long getNewNodeID() { return nodes.size(); };
    int InitMesh(Mesh *mesh);
};

#endif
