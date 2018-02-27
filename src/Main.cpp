#include "Mesh.hpp"
#include "Solver.hpp"
#include "Initializer.hpp"
#include "MeshUtils.hpp"
#include "OutputUtils.hpp"

int main(int argc, char **argv) {
    Mesh *mesh = new Mesh();
    mesh->InitMesh(mesh);
    MeshUtils::calculateNodeNormals(*mesh);
    MeshUtils::calculateVectorsFromCenterToEdges(*mesh);
    Initializer::initialize(*mesh);

    OutputUtils::OutputParaview(mesh, 0);

    Solver solver(mesh);

    for (int i = 0; i < 100; i++) {
        solver.processPhase1();
        solver.processPhase2();
        solver.processPhase3();
        solver.prepareNextStep();
        OutputUtils::OutputParaview(mesh, i + 1);
    }

    return 0;
}
