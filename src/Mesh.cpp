#include "Mesh.hpp"
#include "Parser.hpp"

#include <iostream>

int Mesh::InitMesh(Mesh *mesh) {
    Parser parser;
    parser.LoadNodes(mesh, "triangle/Mesh.1.node");
    parser.LoadEdges(mesh, "triangle/Mesh.1.edge");
    parser.LoadCells(mesh, "triangle/Mesh.1.ele");

    return 0;
}
