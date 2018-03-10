#include "Node.hpp"
#include "Mesh.hpp"

#include <iostream>

Node::Node(Mesh &mesh, double x, double y) {
    data.coords = Vector(x,y);
    ID = mesh.getNewNodeID();
}
