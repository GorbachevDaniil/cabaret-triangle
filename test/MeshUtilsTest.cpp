#include "gtest/gtest.h"

#include "Vector.hpp"
#include "Data.hpp"
#include "Node.hpp"
#include "Edge.hpp"
#include "Cell.hpp"
#include "Mesh.hpp"
#include "MeshUtils.hpp"

#include <vector>

TEST(calculateNodeNormals, Positive) {
    double x1 = 1;
    double y1 = 0;
    double x2 = 3;
    double y2 = 2;
    double length = Vector::length(x1 - x2, y1 - y2);

    Node node1;
    node1.data.coords = Vector(x1, y1);

    Node node2;
    node2.data.coords = Vector(x2, y2);

    Node edgeCenter;

    Edge edge;
    edge.nodeIDs.push_back(0); // node1
    edge.nodeIDs.push_back(1); // node2
    edge.centerNodeID = 2; // edgeCenter
    edge.length = length;

    Cell cell;
    cell.nodeIDs.push_back(0); // node1
    cell.nodeIDs.push_back(1); // node2
    cell.edgeIDs.push_back(0); // edge

    Mesh mesh;
    mesh.nodes.push_back(node1);
    mesh.nodes.push_back(node2);
    mesh.nodes.push_back(edgeCenter);
    mesh.edges.push_back(edge);
    mesh.cells.push_back(cell);

    MeshUtils::calculateNodeNormals(mesh);

    EXPECT_EQ((y2 - y1) / length, mesh.nodes[2].normal.x);
    EXPECT_EQ(-(x2 - x1) / length, mesh.nodes[2].normal.y);
    EXPECT_EQ(1, mesh.cells[0].edgeToNormalDir[0]);
}

TEST(calculateVectorsFromCenterToEdges, Positive) {
    double x1 = 1;
    double y1 = 0;
    double x2 = 3;
    double y2 = 2;
    double length = Vector::length(x1 - x2, y1 - y2);

    Node edgeCenter;
    edgeCenter.data.coords = Vector(x1, y1);

    Node cellCenter;
    cellCenter.data.coords = Vector(x2, y2);

    Edge edge;
    edge.centerNodeID = 0; // edgeCenter

    Cell cell;
    cell.centerNodeID = 1; // cellCenter
    cell.edgeIDs.push_back(0); // edge

    Mesh mesh;
    mesh.nodes.push_back(edgeCenter);
    mesh.nodes.push_back(cellCenter);
    mesh.edges.push_back(edge);
    mesh.cells.push_back(cell);

    MeshUtils::calculateVectorsFromCenterToEdges(mesh);

    EXPECT_EQ((x1 - x2) / length, mesh.cells[0].edgeToTransportDir[0].x);
    EXPECT_EQ((y1 - y2) / length, mesh.cells[0].edgeToTransportDir[0].y);
}
