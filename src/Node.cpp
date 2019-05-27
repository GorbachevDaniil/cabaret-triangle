#include "Node.hpp"
#include "Mesh.hpp"

#include <iostream>

Node::Node(Mesh &mesh, double x, double y, bool used, bool boundNode, 
           bool cellCenterNode,bool isApex) {
    data.coords = Vector(x, y);
    ID = mesh.getNewNodeID();
    this->used = used;
    this->boundNode = boundNode;
    this->cellCenterNode = cellCenterNode;
    this->isApex = isApex;
}
