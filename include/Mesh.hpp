#ifndef Mesh_hpp
#define Mesh_hpp

#include "Cell.hpp"
#include "Edge.hpp"
#include "Node.hpp"

#include <vector>

class Mesh {
public:
  Mesh(int edgeInnerNodesNumber, bool edgeOuterNodesUsed)
      : edgeInnerNodesNumber(edgeInnerNodesNumber),
        edgeOuterNodesUsed(edgeOuterNodesUsed){};

  int edgeInnerNodesNumber;
  bool edgeOuterNodesUsed;

  std::vector<Node> nodes;
  std::vector<Edge> edges;
  std::vector<Cell> cells;
  std::map<std::pair<int, int>, int> mapNodesWithEdge;

  inline long getNewNodeID() { return nodes.size(); };
  int InitMesh(Mesh *mesh);
};

#endif
