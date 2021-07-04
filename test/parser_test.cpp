#include "gtest/gtest.h"

#include "grid/mesh.hpp"
#include "parser.hpp"

class ParserTest : public ::testing::Test
{
public:
    ParserTest() : mesh(1, false) {
        Parser parser;
        parser.LoadNodes(&mesh, "test/resources/triangle/test.node");
        parser.LoadEdges(&mesh, "test/resources/triangle/test.edge");
        parser.LoadCells(&mesh, "test/resources/triangle/test.ele");
    }

protected:
    Mesh mesh;
};

TEST_F(ParserTest, parse_nodes) {
    EXPECT_EQ(-1, mesh.nodes[0].coords.x);
    EXPECT_EQ(-1, mesh.nodes[0].coords.y);
    EXPECT_EQ(-1, mesh.nodes[1].coords.x);
    EXPECT_EQ(1, mesh.nodes[1].coords.y);
    EXPECT_EQ(1, mesh.nodes[2].coords.x);
    EXPECT_EQ(1, mesh.nodes[2].coords.y);
    EXPECT_EQ(1, mesh.nodes[3].coords.x);
    EXPECT_EQ(-1, mesh.nodes[3].coords.y);
}

TEST_F(ParserTest, parse_edges) {
    /* Nodes for every Edge */
    EXPECT_EQ(1, mesh.edges[0].node_ids[0]);
    EXPECT_EQ(4, mesh.edges[0].node_ids[1]);
    EXPECT_EQ(0, mesh.edges[0].node_ids[2]);

    EXPECT_EQ(0, mesh.edges[1].node_ids[0]);
    EXPECT_EQ(5, mesh.edges[1].node_ids[1]);
    EXPECT_EQ(3, mesh.edges[1].node_ids[2]);

    EXPECT_EQ(3, mesh.edges[2].node_ids[0]);
    EXPECT_EQ(6, mesh.edges[2].node_ids[1]);
    EXPECT_EQ(1, mesh.edges[2].node_ids[2]);

    EXPECT_EQ(3, mesh.edges[3].node_ids[0]);
    EXPECT_EQ(7, mesh.edges[3].node_ids[1]);
    EXPECT_EQ(2, mesh.edges[3].node_ids[2]);

    EXPECT_EQ(2, mesh.edges[4].node_ids[0]);
    EXPECT_EQ(8, mesh.edges[4].node_ids[1]);
    EXPECT_EQ(1, mesh.edges[4].node_ids[2]);

    /* Coords for center Nodes */
    EXPECT_EQ((unsigned long) 1, mesh.edges[0].used_node_ids.size());
    long centerID = mesh.edges[0].used_node_ids[0];
    EXPECT_EQ(-1, mesh.nodes[centerID].coords.x);
    EXPECT_EQ(0, mesh.nodes[centerID].coords.y);

    EXPECT_EQ((unsigned long) 1, mesh.edges[1].used_node_ids.size());
    centerID = mesh.edges[1].used_node_ids[0];
    EXPECT_EQ(0, mesh.nodes[centerID].coords.x);
    EXPECT_EQ(-1, mesh.nodes[centerID].coords.y);

    EXPECT_EQ((unsigned long) 1, mesh.edges[2].used_node_ids.size());
    centerID = mesh.edges[2].used_node_ids[0];
    EXPECT_EQ(0, mesh.nodes[centerID].coords.x);
    EXPECT_EQ(0, mesh.nodes[centerID].coords.y);

    EXPECT_EQ((unsigned long) 1, mesh.edges[3].used_node_ids.size());
    centerID = mesh.edges[3].used_node_ids[0];
    EXPECT_EQ(1, mesh.nodes[centerID].coords.x);
    EXPECT_EQ(0, mesh.nodes[centerID].coords.y);

    EXPECT_EQ((unsigned long) 1, mesh.edges[4].used_node_ids.size());
    centerID = mesh.edges[4].used_node_ids[0];
    EXPECT_EQ(0, mesh.nodes[centerID].coords.x);
    EXPECT_EQ(1, mesh.nodes[centerID].coords.y);

    /* Edge's Length */
    EXPECT_EQ(2, mesh.edges[0].length);
    EXPECT_EQ(2, mesh.edges[1].length);
    EXPECT_EQ(sqrt((double)8), mesh.edges[2].length);
    EXPECT_EQ(2, mesh.edges[3].length);
    EXPECT_EQ(2, mesh.edges[4].length);
}

TEST_F(ParserTest, parse_cells) {
    /* Nodes for every Cell */
    EXPECT_EQ(1, mesh.cells[0].node_ids[0]);
    EXPECT_EQ(0, mesh.cells[0].node_ids[1]);
    EXPECT_EQ(3, mesh.cells[0].node_ids[2]);
    EXPECT_EQ(3, mesh.cells[1].node_ids[0]);
    EXPECT_EQ(2, mesh.cells[1].node_ids[1]);
    EXPECT_EQ(1, mesh.cells[1].node_ids[2]);

    /* Coords for center Nodes */
    long centerID = mesh.cells[0].center_node_id;
    EXPECT_EQ(-0.34, floor(mesh.nodes[centerID].coords.x * 100) / 100);
    EXPECT_EQ(-0.34, floor(mesh.nodes[centerID].coords.y * 100) / 100);
    centerID = mesh.cells[1].center_node_id;
    EXPECT_EQ(0.33, floor(mesh.nodes[centerID].coords.x * 100) / 100);
    EXPECT_EQ(0.33, floor(mesh.nodes[centerID].coords.y * 100) / 100);

    /* Edges for every Cell */
    EXPECT_EQ(0, mesh.cells[0].edge_ids[0]);
    EXPECT_EQ(1, mesh.cells[0].edge_ids[1]);
    EXPECT_EQ(2, mesh.cells[0].edge_ids[2]);
    EXPECT_EQ(3, mesh.cells[1].edge_ids[0]);
    EXPECT_EQ(4, mesh.cells[1].edge_ids[1]);
    EXPECT_EQ(2, mesh.cells[1].edge_ids[2]);

    /* Cells for every Edge*/
    EXPECT_EQ(0, mesh.edges[0].cell_ids[0]);
    EXPECT_EQ(0, mesh.edges[1].cell_ids[0]);
    EXPECT_EQ(0, mesh.edges[2].cell_ids[0]);
    EXPECT_EQ(1, mesh.edges[2].cell_ids[1]);
    EXPECT_EQ(1, mesh.edges[3].cell_ids[0]);
    EXPECT_EQ(1, mesh.edges[4].cell_ids[0]);

    /* Cell's Volume */
    EXPECT_EQ(2, mesh.cells[0].volume);
    EXPECT_EQ(2, mesh.cells[1].volume);
}
