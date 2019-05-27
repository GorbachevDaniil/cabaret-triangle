#include "Mesh.hpp"

#include "Parser.hpp"

#include <string>

int Mesh::InitMesh(Mesh *mesh) {
    Parser parser;
    parser.LoadNodes(mesh, "mesh_generation/triangle/Mesh.1.node");
    parser.LoadEdges(mesh, "mesh_generation/triangle/Mesh.1.edge");
    parser.LoadCells(mesh, "mesh_generation/triangle/Mesh.1.ele");

    return 0;
}