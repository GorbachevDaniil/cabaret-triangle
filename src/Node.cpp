#include "Node.hpp"
#include "Mesh.hpp"

#include <iostream>

Node::Node(Mesh &mesh, double x, double y, bool used) {
    data.coords = Vector(x,y);
    ID = mesh.getNewNodeID();
    this->used = used;
}
