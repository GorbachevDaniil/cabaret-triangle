#include "Mesh.hpp"

#include "Parser.hpp"

#include <string>

int Mesh::InitMesh(Mesh *mesh) {
    Parser parser;
    parser.LoadNodes(mesh, "triangle_regular/Mesh.node");
    parser.LoadEdges(mesh, "triangle_regular/Mesh.edge");
    parser.LoadCells(mesh, "triangle_regular/Mesh.ele");

    return 0;
}
