#include "gtest/gtest.h"

#include "Node.hpp"
#include "Edge.hpp"
#include "Cell.hpp"
#include "Mesh.hpp"

#include <vector>
#include <iostream>

TEST(getNextNodeIDTest, ReturnTheNextValueInACircularManner) {
    Cell cell;
    cell.nodeIDs.push_back(10);
    cell.nodeIDs.push_back(11);
    cell.nodeIDs.push_back(12);

    EXPECT_EQ(11, cell.getNextNodeID(0));
    EXPECT_EQ(12, cell.getNextNodeID(1));
    EXPECT_EQ(10, cell.getNextNodeID(2));
}

TEST(getNextNodeIDTest, DieIfIncomingNodeIDPosMoreThanSize) {
    Cell cell;

    ASSERT_DEATH(cell.getNextNodeID(1), "");
}

TEST(getPrevNodeIDTest, ReturnTheNextValueInACircularManner) {
    Cell cell;
    cell.nodeIDs.push_back(10);
    cell.nodeIDs.push_back(11);
    cell.nodeIDs.push_back(12);

    EXPECT_EQ(12, cell.getPrevNodeID(0));
    EXPECT_EQ(10, cell.getPrevNodeID(1));
    EXPECT_EQ(11, cell.getPrevNodeID(2));
}

TEST(getPrevNodeIDTest, DieIfIncomingNodeIDPosMoreThanSize) {
    Cell cell;

    ASSERT_DEATH(cell.getNextNodeID(1), "");
}

TEST(getEdgeOrderedNodeIDs, UnorderedNodeIDs) {
    Cell cell;
    cell.nodeIDs.push_back(10);
    cell.nodeIDs.push_back(11);
    cell.nodeIDs.push_back(12);

    std::vector<long> unorderedNodeIDs;
    unorderedNodeIDs.push_back(10);
    unorderedNodeIDs.push_back(12);

    std::vector<long> orderedNodeIDs = cell.getEdgeOrderedNodeIDs(unorderedNodeIDs);
    EXPECT_EQ(12, orderedNodeIDs[0]);
    EXPECT_EQ(10, orderedNodeIDs[1]);
}

TEST(getEdgeOrderedNodeIDs, OrderedNodeIDs) {
    Cell cell;
    cell.nodeIDs.push_back(10);
    cell.nodeIDs.push_back(11);
    cell.nodeIDs.push_back(12);

    std::vector<long> unorderedNodeIDs;
    unorderedNodeIDs.push_back(12);
    unorderedNodeIDs.push_back(10);

    std::vector<long> orderedNodeIDs = cell.getEdgeOrderedNodeIDs(unorderedNodeIDs);
    EXPECT_EQ(12, orderedNodeIDs[0]);
    EXPECT_EQ(10, orderedNodeIDs[1]);
}

TEST(getEdgeOrderedNodeIDs, SizeMustBeTwoForIncomingVector) {
    Cell cell;

    std::vector<long> unorderedNodeIDs;
    unorderedNodeIDs.push_back(10);
    unorderedNodeIDs.push_back(11);
    unorderedNodeIDs.push_back(12);

    ASSERT_DEATH(cell.getEdgeOrderedNodeIDs(unorderedNodeIDs), "");

    unorderedNodeIDs.clear();
    unorderedNodeIDs.push_back(10);

    ASSERT_DEATH(cell.getEdgeOrderedNodeIDs(unorderedNodeIDs), "");
}

TEST(getEdgeOrderedNodeIDs, SizeMustBeTwoForOutcomingVector) {
    Cell cell;
    cell.nodeIDs.push_back(10);
    cell.nodeIDs.push_back(11);
    cell.nodeIDs.push_back(12);

    std::vector<long> unorderedNodeIDs;
    unorderedNodeIDs.push_back(11);
    unorderedNodeIDs.push_back(11);

    ASSERT_DEATH(cell.getEdgeOrderedNodeIDs(unorderedNodeIDs), "");
}

TEST(cellCreationInnerNodes, OneInnerNodeDontHaveOpposite) {
    double x1 = 0;
    double y1 = 0;

    double x2 = 1;
    double y2 = 0;

    double x3 = 0;
    double y3 = 1;

    Node node1;
    node1.ID = 0;
    node1.used = false;
    node1.data.coords = Vector(x1, y1);

    Node node2;
    node2.ID = 1;
    node2.used = false;
    node2.data.coords = Vector(x2, y2);

    Node node3;
    node3.ID = 2;
    node3.used = false;
    node3.data.coords = Vector(x3, y3);

    Mesh mesh;
    mesh.nodes.push_back(node1);
    mesh.nodes.push_back(node2);
    mesh.nodes.push_back(node3);

    Edge edge1 = Edge(mesh, 0, node1.ID, node2.ID, false, 1);
    Edge edge2 = Edge(mesh, 1, node2.ID, node3.ID, false, 1);
    Edge edge3 = Edge(mesh, 2, node3.ID, node1.ID, false, 1);

    mesh.edges.push_back(edge1);
    mesh.edges.push_back(edge2);
    mesh.edges.push_back(edge3);

    Cell cell = Cell(mesh, 0, node1.ID, node2.ID, node3.ID);

    EXPECT_EQ((unsigned long) 7, mesh.nodes.size());

    Node *edge1Node = &mesh.nodes[3];
    EXPECT_EQ((x1 + x2) / 2, edge1Node->data.coords.x);
    EXPECT_EQ((y1 + y2) / 2, edge1Node->data.coords.y);
    EXPECT_EQ(true, edge1Node->used);
    EXPECT_EQ(-1, cell.nodeIDToOppositeNodeID[edge1Node->ID]);

    Node *edge2Node = &mesh.nodes[4];
    EXPECT_EQ((x2 + x3) / 2, edge2Node->data.coords.x);
    EXPECT_EQ((y2 + y3) / 2, edge2Node->data.coords.y);
    EXPECT_EQ(true, edge2Node->used);
    EXPECT_EQ(-1, cell.nodeIDToOppositeNodeID[edge2Node->ID]);

    Node *edge3Node = &mesh.nodes[5];
    EXPECT_EQ((x1 + x3) / 2, edge3Node->data.coords.x);
    EXPECT_EQ((y1 + y3) / 2, edge3Node->data.coords.y);
    EXPECT_EQ(true, edge3Node->used);
    EXPECT_EQ(-1, cell.nodeIDToOppositeNodeID[edge3Node->ID]);
}

TEST(cellCreationInnerNodes, TwoInnerNodesHaveOppositeWithPositiveCoordsLeftTriangle) {
    double x1 = 0;
    double y1 = 0;

    double x2 = 1;
    double y2 = 0;

    double x3 = 0;
    double y3 = 1;

    Node node1;
    node1.ID = 0;
    node1.used = false;
    node1.data.coords = Vector(x1, y1);

    Node node2;
    node2.ID = 1;
    node2.used = false;
    node2.data.coords = Vector(x2, y2);

    Node node3;
    node3.ID = 2;
    node3.used = false;
    node3.data.coords = Vector(x3, y3);

    Mesh mesh;
    mesh.nodes.push_back(node1);
    mesh.nodes.push_back(node2);
    mesh.nodes.push_back(node3);

    Edge edge1 = Edge(mesh, 0, node1.ID, node2.ID, false, 2);
    Edge edge2 = Edge(mesh, 1, node2.ID, node3.ID, false, 2);
    Edge edge3 = Edge(mesh, 2, node3.ID, node1.ID, false, 2);

    mesh.edges.push_back(edge1);
    mesh.edges.push_back(edge2);
    mesh.edges.push_back(edge3);

    Cell cell = Cell(mesh, 0, node1.ID, node2.ID, node3.ID);

    EXPECT_EQ((unsigned long) 10, mesh.nodes.size());

    Node *edge1Node1 = &mesh.nodes[3];
    EXPECT_DOUBLE_EQ(1. / 3., edge1Node1->data.coords.x);
    EXPECT_DOUBLE_EQ(0, edge1Node1->data.coords.y);
    EXPECT_EQ(true, edge1Node1->used);
    EXPECT_EQ(6, cell.nodeIDToOppositeNodeID[edge1Node1->ID]);

    Node *edge1Node2 = &mesh.nodes[4];
    EXPECT_DOUBLE_EQ(2. / 3., edge1Node2->data.coords.x);
    EXPECT_DOUBLE_EQ(0, edge1Node2->data.coords.y);
    EXPECT_EQ(true, edge1Node2->used);
    EXPECT_EQ(7, cell.nodeIDToOppositeNodeID[edge1Node2->ID]);

    Node *edge2Node1 = &mesh.nodes[5];
    EXPECT_DOUBLE_EQ(2. / 3., edge2Node1->data.coords.x);
    EXPECT_DOUBLE_EQ(1. / 3., edge2Node1->data.coords.y);
    EXPECT_EQ(true, edge2Node1->used);
    EXPECT_EQ(8, cell.nodeIDToOppositeNodeID[edge2Node1->ID]);

    Node *edge2Node2 = &mesh.nodes[6];
    EXPECT_DOUBLE_EQ(1. / 3., edge2Node2->data.coords.x);
    EXPECT_DOUBLE_EQ(2. / 3., edge2Node2->data.coords.y);
    EXPECT_EQ(true, edge2Node2->used);
    EXPECT_EQ(3, cell.nodeIDToOppositeNodeID[edge2Node2->ID]);

    Node *edge3Node1 = &mesh.nodes[7];
    EXPECT_DOUBLE_EQ(0, edge3Node1->data.coords.x);
    EXPECT_DOUBLE_EQ(2. / 3., edge3Node1->data.coords.y);
    EXPECT_EQ(true, edge3Node1->used);
    EXPECT_EQ(4, cell.nodeIDToOppositeNodeID[edge3Node1->ID]);

    Node *edge3Node2 = &mesh.nodes[8];
    EXPECT_DOUBLE_EQ(0, edge3Node2->data.coords.x);
    EXPECT_DOUBLE_EQ(1. / 3., edge3Node2->data.coords.y);
    EXPECT_EQ(true, edge3Node2->used);
    EXPECT_EQ(5, cell.nodeIDToOppositeNodeID[edge3Node2->ID]);
}

TEST(cellCreationInnerNodes, TwoInnerNodesHaveOppositeWithPositiveCoordsRightTriangle) {
    double x1 = 1;
    double y1 = 0;

    double x2 = 1;
    double y2 = 1;

    double x3 = 0;
    double y3 = 1;

    Node node1;
    node1.ID = 0;
    node1.used = false;
    node1.data.coords = Vector(x1, y1);

    Node node2;
    node2.ID = 1;
    node2.used = false;
    node2.data.coords = Vector(x2, y2);

    Node node3;
    node3.ID = 2;
    node3.used = false;
    node3.data.coords = Vector(x3, y3);

    Mesh mesh;
    mesh.nodes.push_back(node1);
    mesh.nodes.push_back(node2);
    mesh.nodes.push_back(node3);

    Edge edge1 = Edge(mesh, 0, node1.ID, node2.ID, false, 2);
    Edge edge2 = Edge(mesh, 1, node2.ID, node3.ID, false, 2);
    Edge edge3 = Edge(mesh, 2, node3.ID, node1.ID, false, 2);

    mesh.edges.push_back(edge1);
    mesh.edges.push_back(edge2);
    mesh.edges.push_back(edge3);

    Cell cell = Cell(mesh, 0, node1.ID, node2.ID, node3.ID);

    EXPECT_EQ((unsigned long) 10, mesh.nodes.size());

    Node *edge1Node1 = &mesh.nodes[3];
    EXPECT_DOUBLE_EQ(1, edge1Node1->data.coords.x);
    EXPECT_DOUBLE_EQ(1. / 3., edge1Node1->data.coords.y);
    EXPECT_EQ(true, edge1Node1->used);
    EXPECT_EQ(6, cell.nodeIDToOppositeNodeID[edge1Node1->ID]);

    Node *edge1Node2 = &mesh.nodes[4];
    EXPECT_DOUBLE_EQ(1, edge1Node2->data.coords.x);
    EXPECT_DOUBLE_EQ(2. / 3., edge1Node2->data.coords.y);
    EXPECT_EQ(true, edge1Node2->used);
    EXPECT_EQ(7, cell.nodeIDToOppositeNodeID[edge1Node2->ID]);

    Node *edge2Node1 = &mesh.nodes[5];
    EXPECT_DOUBLE_EQ(2. / 3., edge2Node1->data.coords.x);
    EXPECT_DOUBLE_EQ(1, edge2Node1->data.coords.y);
    EXPECT_EQ(true, edge2Node1->used);
    EXPECT_EQ(8, cell.nodeIDToOppositeNodeID[edge2Node1->ID]);

    Node *edge2Node2 = &mesh.nodes[6];
    EXPECT_DOUBLE_EQ(1. / 3., edge2Node2->data.coords.x);
    EXPECT_DOUBLE_EQ(1, edge2Node2->data.coords.y);
    EXPECT_EQ(true, edge2Node2->used);
    EXPECT_EQ(3, cell.nodeIDToOppositeNodeID[edge2Node2->ID]);

    Node *edge3Node1 = &mesh.nodes[7];
    EXPECT_DOUBLE_EQ(1. / 3., edge3Node1->data.coords.x);
    EXPECT_DOUBLE_EQ(2. / 3., edge3Node1->data.coords.y);
    EXPECT_EQ(true, edge3Node1->used);
    EXPECT_EQ(4, cell.nodeIDToOppositeNodeID[edge3Node1->ID]);

    Node *edge3Node2 = &mesh.nodes[8];
    EXPECT_DOUBLE_EQ(2. / 3., edge3Node2->data.coords.x);
    EXPECT_DOUBLE_EQ(1. / 3., edge3Node2->data.coords.y);
    EXPECT_EQ(true, edge3Node2->used);
    EXPECT_EQ(5, cell.nodeIDToOppositeNodeID[edge3Node2->ID]);
}

TEST(cellCreationInnerNodes, TwoInnerNodesHaveOppositeWithNegativeCoords) {
    double x1 = 0;
    double y1 = -1;

    double x2 = 0;
    double y2 = 0;

    double x3 = -1;
    double y3 = 0;

    Node node1;
    node1.ID = 0;
    node1.used = false;
    node1.data.coords = Vector(x1, y1);

    Node node2;
    node2.ID = 1;
    node2.used = false;
    node2.data.coords = Vector(x2, y2);

    Node node3;
    node3.ID = 2;
    node3.used = false;
    node3.data.coords = Vector(x3, y3);

    Mesh mesh;
    mesh.nodes.push_back(node1);
    mesh.nodes.push_back(node2);
    mesh.nodes.push_back(node3);

    Edge edge1 = Edge(mesh, 0, node1.ID, node2.ID, false, 2);
    Edge edge2 = Edge(mesh, 1, node2.ID, node3.ID, false, 2);
    Edge edge3 = Edge(mesh, 2, node3.ID, node1.ID, false, 2);

    mesh.edges.push_back(edge1);
    mesh.edges.push_back(edge2);
    mesh.edges.push_back(edge3);

    Cell cell = Cell(mesh, 0, node1.ID, node2.ID, node3.ID);

    EXPECT_EQ((unsigned long) 10, mesh.nodes.size());

    Node *edge1Node1 = &mesh.nodes[3];
    EXPECT_DOUBLE_EQ(0, edge1Node1->data.coords.x);
    EXPECT_DOUBLE_EQ(-2. / 3., edge1Node1->data.coords.y);
    EXPECT_EQ(true, edge1Node1->used);
    EXPECT_EQ(6, cell.nodeIDToOppositeNodeID[edge1Node1->ID]);

    Node *edge1Node2 = &mesh.nodes[4];
    EXPECT_DOUBLE_EQ(0, edge1Node2->data.coords.x);
    EXPECT_DOUBLE_EQ(-1. / 3., edge1Node2->data.coords.y);
    EXPECT_EQ(true, edge1Node2->used);
    EXPECT_EQ(7, cell.nodeIDToOppositeNodeID[edge1Node2->ID]);

    Node *edge2Node1 = &mesh.nodes[5];
    EXPECT_DOUBLE_EQ(-1. / 3., edge2Node1->data.coords.x);
    EXPECT_DOUBLE_EQ(0, edge2Node1->data.coords.y);
    EXPECT_EQ(true, edge2Node1->used);
    EXPECT_EQ(8, cell.nodeIDToOppositeNodeID[edge2Node1->ID]);

    Node *edge2Node2 = &mesh.nodes[6];
    EXPECT_DOUBLE_EQ(-2. / 3., edge2Node2->data.coords.x);
    EXPECT_DOUBLE_EQ(0, edge2Node2->data.coords.y);
    EXPECT_EQ(true, edge2Node2->used);
    EXPECT_EQ(3, cell.nodeIDToOppositeNodeID[edge2Node2->ID]);

    Node *edge3Node1 = &mesh.nodes[7];
    EXPECT_DOUBLE_EQ(-2. / 3., edge3Node1->data.coords.x);
    EXPECT_DOUBLE_EQ(-1. / 3., edge3Node1->data.coords.y);
    EXPECT_EQ(true, edge3Node1->used);
    EXPECT_EQ(4, cell.nodeIDToOppositeNodeID[edge3Node1->ID]);

    Node *edge3Node2 = &mesh.nodes[8];
    EXPECT_DOUBLE_EQ(-1. / 3., edge3Node2->data.coords.x);
    EXPECT_DOUBLE_EQ(-2. / 3., edge3Node2->data.coords.y);
    EXPECT_EQ(true, edge3Node2->used);
    EXPECT_EQ(5, cell.nodeIDToOppositeNodeID[edge3Node2->ID]);
}

TEST(cellCreationInnerNodes, TwoInnerNodesHaveOppositeWithWrongOrderOnEdge) {
    double x1 = 1;
    double y1 = 0;

    double x2 = 1;
    double y2 = 1;

    double x3 = 0;
    double y3 = 1;

    Node node1;
    node1.ID = 0;
    node1.used = false;
    node1.data.coords = Vector(x1, y1);

    Node node2;
    node2.ID = 1;
    node2.used = false;
    node2.data.coords = Vector(x2, y2);

    Node node3;
    node3.ID = 2;
    node3.used = false;
    node3.data.coords = Vector(x3, y3);

    Mesh mesh;
    mesh.nodes.push_back(node1);
    mesh.nodes.push_back(node2);
    mesh.nodes.push_back(node3);

    Edge edge1 = Edge(mesh, 0, node1.ID, node2.ID, false, 2);
    Edge edge2 = Edge(mesh, 1, node3.ID, node2.ID, false, 2);
    Edge edge3 = Edge(mesh, 2, node3.ID, node1.ID, false, 2);

    mesh.edges.push_back(edge1);
    mesh.edges.push_back(edge2);
    mesh.edges.push_back(edge3);

    Cell cell = Cell(mesh, 0, node1.ID, node2.ID, node3.ID);

    EXPECT_EQ((unsigned long) 10, mesh.nodes.size());

    Node *edge1Node1 = &mesh.nodes[3];
    EXPECT_DOUBLE_EQ(1, edge1Node1->data.coords.x);
    EXPECT_DOUBLE_EQ(1. / 3., edge1Node1->data.coords.y);
    EXPECT_EQ(true, edge1Node1->used);
    EXPECT_EQ(5, cell.nodeIDToOppositeNodeID[edge1Node1->ID]);

    Node *edge1Node2 = &mesh.nodes[4];
    EXPECT_DOUBLE_EQ(1, edge1Node2->data.coords.x);
    EXPECT_DOUBLE_EQ(2. / 3., edge1Node2->data.coords.y);
    EXPECT_EQ(true, edge1Node2->used);
    EXPECT_EQ(7, cell.nodeIDToOppositeNodeID[edge1Node2->ID]);

    Node *edge2Node1 = &mesh.nodes[5];
    EXPECT_DOUBLE_EQ(1. / 3., edge2Node1->data.coords.x);
    EXPECT_DOUBLE_EQ(1, edge2Node1->data.coords.y);
    EXPECT_EQ(true, edge2Node1->used);
    EXPECT_EQ(3, cell.nodeIDToOppositeNodeID[edge2Node1->ID]);

    Node *edge2Node2 = &mesh.nodes[6];
    EXPECT_DOUBLE_EQ(2. / 3., edge2Node2->data.coords.x);
    EXPECT_DOUBLE_EQ(1, edge2Node2->data.coords.y);
    EXPECT_EQ(true, edge2Node2->used);
    EXPECT_EQ(8, cell.nodeIDToOppositeNodeID[edge2Node2->ID]);

    Node *edge3Node1 = &mesh.nodes[7];
    EXPECT_DOUBLE_EQ(1. / 3., edge3Node1->data.coords.x);
    EXPECT_DOUBLE_EQ(2. / 3., edge3Node1->data.coords.y);
    EXPECT_EQ(true, edge3Node1->used);
    EXPECT_EQ(4, cell.nodeIDToOppositeNodeID[edge3Node1->ID]);

    Node *edge3Node2 = &mesh.nodes[8];
    EXPECT_DOUBLE_EQ(2. / 3., edge3Node2->data.coords.x);
    EXPECT_DOUBLE_EQ(1. / 3., edge3Node2->data.coords.y);
    EXPECT_EQ(true, edge3Node2->used);
    EXPECT_EQ(6, cell.nodeIDToOppositeNodeID[edge3Node2->ID]);
}