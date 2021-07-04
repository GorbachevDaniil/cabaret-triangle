#include "grid/node.hpp"

#include <iostream>

#include "grid/mesh.hpp"

Node::Node(Mesh &mesh, double x, double y, bool used, bool is_bound) {
    coords = Vector(x, y);
    id = mesh.get_new_node_id();
    this->used = used;
    this->is_bound = is_bound;
}
