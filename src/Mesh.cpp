#include "Mesh.hpp"

#include "Parser.hpp"

#include <string>

int Mesh::InitMesh(Mesh *mesh) {
    Parser parser;
    parser.LoadNodes(mesh, "mesh/chaos/square 50x50/mesh.1.node");
    parser.LoadEdges(mesh, "mesh/chaos/square 50x50/mesh.1.edge");
    parser.LoadCells(mesh, "mesh/chaos/square 50x50/mesh.1.ele");
    
    long nodesSize = mesh->nodes.size();

    mesh->s0 = std::vector<std::vector<double>>(nodesSize);
    mesh->s1 = std::vector<std::vector<double>>(nodesSize);
    mesh->s2 = std::vector<std::vector<double>>(nodesSize);

    mesh->v0 = std::vector<std::vector<Vector>>(nodesSize);
    mesh->v1 = std::vector<std::vector<Vector>>(nodesSize);
    mesh->v2 = std::vector<std::vector<Vector>>(nodesSize);

    return 0;
}