#include "Mesh.hpp"

#include "Parser.hpp"

#include <string>

int Mesh::InitMesh(Mesh *mesh) {
    Parser parser;
    parser.LoadNodes(mesh, "mesh_generation/triangle_regular/Mesh.node");
    parser.LoadEdges(mesh, "mesh_generation/triangle_regular/Mesh.edge");
    parser.LoadCells(mesh, "mesh_generation/triangle_regular/Mesh.ele");

    return 0;
}
