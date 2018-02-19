#include "Node.hpp"
#include "Mesh.hpp"

#include <iostream>

Node::Node(Mesh &mesh, double x, double y){
	Vector vec;
	vec.set(x,y);
	data.coords = vec;
	ID = mesh.getNewNodeID();
}
