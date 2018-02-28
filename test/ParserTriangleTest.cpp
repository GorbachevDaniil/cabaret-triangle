#include "gtest/gtest.h"

#include "Parser.hpp"
#include "Mesh.hpp"


TEST(ParserTriangleTest, LoadNodes) {
    Parser parser;
    Mesh mesh;
    parser.LoadNodes(&mesh, "test/resources/triangle/Mesh.1.node");

    EXPECT_EQ( -1, mesh.nodes[0].data.coords.x);
    EXPECT_EQ( -1, mesh.nodes[0].data.coords.y);
    EXPECT_EQ( -1, mesh.nodes[1].data.coords.x);
    EXPECT_EQ( 1, mesh.nodes[1].data.coords.y);
    EXPECT_EQ( 1, mesh.nodes[2].data.coords.x);
    EXPECT_EQ( 1, mesh.nodes[2].data.coords.y);
    EXPECT_EQ( 1, mesh.nodes[3].data.coords.x);
    EXPECT_EQ( -1, mesh.nodes[3].data.coords.y);
}   

TEST(ParserTriangleTest, LoadEdges) {
    Parser parser;
    Mesh mesh;
    parser.LoadNodes(&mesh, "test/resources/triangle/Mesh.1.node");
    parser.LoadEdges(&mesh, "test/resources/triangle/Mesh.1.edge");

    /* Nodes for every Edge */
    EXPECT_EQ( 1, mesh.edges[0].nodeIDs[0]);
    EXPECT_EQ( 0, mesh.edges[0].nodeIDs[1]);
    EXPECT_EQ( 0, mesh.edges[1].nodeIDs[0]);
    EXPECT_EQ( 3, mesh.edges[1].nodeIDs[1]);
    EXPECT_EQ( 3, mesh.edges[2].nodeIDs[0]);
    EXPECT_EQ( 1, mesh.edges[2].nodeIDs[1]);
    EXPECT_EQ( 3, mesh.edges[3].nodeIDs[0]);
    EXPECT_EQ( 2, mesh.edges[3].nodeIDs[1]);
    EXPECT_EQ( 2, mesh.edges[4].nodeIDs[0]);
    EXPECT_EQ( 1, mesh.edges[4].nodeIDs[1]);

    /* Coords for center Nodes */
    long centerID = mesh.edges[0].centerNodeID;
    EXPECT_EQ( -1, mesh.nodes[centerID].data.coords.x);
    EXPECT_EQ( 0, mesh.nodes[centerID].data.coords.y);
    centerID = mesh.edges[1].centerNodeID;
    EXPECT_EQ( 0, mesh.nodes[centerID].data.coords.x);
    EXPECT_EQ( -1, mesh.nodes[centerID].data.coords.y);
    centerID = mesh.edges[2].centerNodeID;
    EXPECT_EQ( 0, mesh.nodes[centerID].data.coords.x);
    EXPECT_EQ( 0, mesh.nodes[centerID].data.coords.y);
    centerID = mesh.edges[3].centerNodeID;
    EXPECT_EQ( 1, mesh.nodes[centerID].data.coords.x);
    EXPECT_EQ( 0, mesh.nodes[centerID].data.coords.y);
    centerID = mesh.edges[4].centerNodeID;
    EXPECT_EQ( 0, mesh.nodes[centerID].data.coords.x);
    EXPECT_EQ( 1, mesh.nodes[centerID].data.coords.y);

    /* Edge's Length */
    EXPECT_EQ( 2, mesh.edges[0].length);
    EXPECT_EQ( 2, mesh.edges[1].length);
    EXPECT_EQ( sqrt(8), mesh.edges[2].length);
    EXPECT_EQ( 2, mesh.edges[3].length);
    EXPECT_EQ( 2, mesh.edges[4].length);
}     


TEST(ParserTriangleTest, LoadCells) {
    Parser parser;
    Mesh mesh;
    parser.LoadNodes(&mesh, "test/resources/triangle/Mesh.1.node");
    parser.LoadEdges(&mesh, "test/resources/triangle/Mesh.1.edge");
    parser.LoadCells(&mesh, "test/resources/triangle/Mesh.1.ele");

    /* Nodes for every Cell */
    EXPECT_EQ( 1, mesh.cells[0].nodeIDs[0]);
    EXPECT_EQ( 0, mesh.cells[0].nodeIDs[1]);
    EXPECT_EQ( 3, mesh.cells[0].nodeIDs[2]);
    EXPECT_EQ( 3, mesh.cells[1].nodeIDs[0]);
    EXPECT_EQ( 2, mesh.cells[1].nodeIDs[1]);
    EXPECT_EQ( 1, mesh.cells[1].nodeIDs[2]);

    /* Coords for center Nodes */
    long centerID = mesh.cells[0].centerNodeID;
    EXPECT_EQ( -0.34, floor(mesh.nodes[centerID].data.coords.x*100)/100);
    EXPECT_EQ( -0.34, floor(mesh.nodes[centerID].data.coords.y*100)/100);
    centerID = mesh.cells[1].centerNodeID;
    EXPECT_EQ( 0.33, floor(mesh.nodes[centerID].data.coords.x*100)/100);
    EXPECT_EQ( 0.33, floor(mesh.nodes[centerID].data.coords.y*100)/100);

    /* Edges for every Cell */
    EXPECT_EQ( 0, mesh.cells[0].edgeIDs[0]);
    EXPECT_EQ( 1, mesh.cells[0].edgeIDs[1]);
    EXPECT_EQ( 2, mesh.cells[0].edgeIDs[2]);
    EXPECT_EQ( 3, mesh.cells[1].edgeIDs[0]);
    EXPECT_EQ( 4, mesh.cells[1].edgeIDs[1]);
    EXPECT_EQ( 2, mesh.cells[1].edgeIDs[2]);

    /* Cells for every Edge*/
    EXPECT_EQ( 0, mesh.edges[0].cellIDs[0]);
    EXPECT_EQ( 0, mesh.edges[1].cellIDs[0]);
    EXPECT_EQ( 0, mesh.edges[2].cellIDs[0]);
    EXPECT_EQ( 1, mesh.edges[2].cellIDs[1]);
    EXPECT_EQ( 1, mesh.edges[3].cellIDs[0]);
    EXPECT_EQ( 1, mesh.edges[4].cellIDs[0]);

    /* Cell's Volume */
    EXPECT_EQ( 2, mesh.cells[0].volume);
    EXPECT_EQ( 2, mesh.cells[1].volume);
}     



