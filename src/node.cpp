#include "node.hpp"

#include <iostream>

#include "mesh.hpp"

Node::Node(Mesh &mesh, double x, double y, bool used, bool boundNode,
           bool cellCenterNode, bool isApex) {
    data.coords = Vector(x, y);
    ID = mesh.getNewNodeID();
    data.nodeID = ID;
    data.mesh = &mesh;
    this->used = used;
    this->boundNode = boundNode;
    this->cellCenterNode = cellCenterNode;
    this->isApex = isApex;
}
