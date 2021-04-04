#include "gtest/gtest.h"

#include <vector>

#include "vector.hpp"
#include "grid/node.hpp"
#include "grid/edge.hpp"
#include "grid/mesh.hpp"

TEST(edgeInnerNodes, EdgeWithOneInnerNodeMustContainThreeNodes) {
    double x1 = 1;
    double y1 = 0;
    double x2 = 3;
    double y2 = 2;

    Node node1;
    node1.ID = 0;
    node1.data.coords = Vector(x1, y1);

    Node node2;
    node2.ID = 1;
    node2.data.coords = Vector(x2, y2);

    Mesh mesh(1, false);
    mesh.nodes.push_back(node1);
    mesh.nodes.push_back(node2);

    Edge edge = Edge(mesh, 0, node1.ID, node2.ID, false, 1);

    EXPECT_EQ((unsigned long) 3, mesh.nodes.size());

    Node *node3 = &mesh.nodes[2];
    ASSERT_DOUBLE_EQ((x1 + x2) / 2, node3->data.coords.x);
    ASSERT_DOUBLE_EQ((y1 + y2) / 2, node3->data.coords.y);
    EXPECT_EQ(true, node3->used);
}

TEST(edgeInnerNodes, EdgeWithTwoInnerNodesMustContainFourNodes) {
    double x1 = 1;
    double y1 = 0;
    double x2 = 3;
    double y2 = 2;

    Node node1;
    node1.ID = 0;
    node1.data.coords = Vector(x1, y1);

    Node node2;
    node2.ID = 1;
    node2.data.coords = Vector(x2, y2);

    Mesh mesh(2, false);
    mesh.nodes.push_back(node1);
    mesh.nodes.push_back(node2);

    Edge edge = Edge(mesh, 0, node1.ID, node2.ID, false, 2);

    EXPECT_EQ((unsigned long) 4, mesh.nodes.size());

    Node *node3 = &mesh.nodes[2];
    EXPECT_DOUBLE_EQ(5. / 3., node3->data.coords.x);
    EXPECT_DOUBLE_EQ(2. / 3., node3->data.coords.y);
    EXPECT_EQ(true, node3->used);

    Node *node4 = &mesh.nodes[3];
    EXPECT_DOUBLE_EQ(7. / 3., node4->data.coords.x);
    EXPECT_DOUBLE_EQ(4. / 3., node4->data.coords.y);
    EXPECT_EQ(true, node4->used);
}