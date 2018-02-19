#include "gtest/gtest.h"

#include "Vector.hpp"
#include "Data.hpp"
#include "Node.hpp"
#include "Edge.hpp"
#include "Cell.hpp"
#include "Mesh.hpp"
#include "MeshUtils.hpp"

#include <vector>
#include <iostream>

TEST(calculateNormals, Positive) {
    double x1 = 1;
    double y1 = 0;
    double x2 = 3;
    double y2 = 2;

    Node node1;
    node1.data.coords.set(x1, y1);

    Node node2;
    node2.data.coords.set(x2, y2);

    Node nodeCenter;

    Edge edge;
    edge.nodeIDs.push_back(0); // node1
    edge.nodeIDs.push_back(1); // node2
    edge.centerNodeID = 2; // nodeCenter
    edge.length = Vector::Length(x1 - x2, y1 - y2);

    Cell cell;
    cell.edgeIDs.push_back(0); // edge
    cell.nodeIDs.push_back(0); // node1
    cell.nodeIDs.push_back(1); // node2

    Mesh mesh;
    mesh.nodes.push_back(node1);
    mesh.nodes.push_back(node2);
    mesh.nodes.push_back(nodeCenter);
    mesh.edges.push_back(edge);
    mesh.cells.push_back(cell);

    MeshUtils::calculateNormals(mesh);

    EXPECT_EQ((y2 - y1) / edge.length, mesh.nodes[2].normal.x); 
    EXPECT_EQ(-(x2 - x1) / edge.length, mesh.nodes[2].normal.y);
}