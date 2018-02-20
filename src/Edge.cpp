#include "Edge.hpp"
#include "Node.hpp"
#include "Mesh.hpp"

#include <iostream>
#include <math.h>

Edge::Edge(Mesh &mesh, long id, long start_node_id, long end_node_id, bool boundary){
	
	double start_node_x, start_node_y, end_node_x, end_node_y;

	start_node_x = mesh.nodes[start_node_id].data.coords.x;
	start_node_y = mesh.nodes[start_node_id].data.coords.y;
	end_node_x = mesh.nodes[end_node_id].data.coords.x;
	end_node_y = mesh.nodes[end_node_id].data.coords.y;

	ID = id;
	nodeIDs.push_back(start_node_id);
	nodeIDs.push_back(end_node_id);
	boundEdge = boundary;
	length = sqrt(pow((end_node_x-start_node_x), 2) + pow((end_node_y-start_node_y), 2));

	Node *node = new Node(mesh, (start_node_x + end_node_x) / 2, (start_node_y + end_node_y) / 2);
	mesh.nodes.push_back(*node);
	centerNodeID = node->ID;
}
