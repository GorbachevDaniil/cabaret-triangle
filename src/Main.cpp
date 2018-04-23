#include "Mesh.hpp"
#include "Solver.hpp"
#include "Initializer.hpp"
#include "MeshUtils.hpp"
#include "OutputUtils.hpp"
#include "Parameters.hpp"

#include <iostream>

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Number of steps should be specified" << std::endl;
        return 1;
    }
    int steps = std::atoi(argv[1]);

    Mesh *mesh = new Mesh();
    mesh->InitMesh(mesh);
    std::cout << "number of cells = " << mesh->cells.size() << std::endl;

    MeshUtils::calculateNodeNormals(*mesh);
    MeshUtils::calculateVectorsFromCenterToEdges(*mesh);
    Initializer::initialize(*mesh);

    OutputUtils::OutputParaview(mesh, 0);

    Solver solver(mesh);
    double tau = solver.calculateTau();
    std::cout << "tau = " << tau << std::endl;

    for (int i = 0; i < steps; i++) {
        solver.processPhase1(tau);
        solver.processPhase2(tau);
        solver.processPhase3(tau);
        solver.prepareNextStep();
        if (steps % Parameters::HARMONIZATION_STEP == 0) {
            solver.harmonise();
        }

        OutputUtils::OutputParaview(mesh, i + 1);
    }

    return 0;
}
