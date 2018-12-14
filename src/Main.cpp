#include "Mesh.hpp"
#include "Solver.hpp"
#include "Initializer.hpp"
#include "MeshUtils.hpp"
#include "OutputUtils.hpp"
#include "Parameters.hpp"

#include <iostream>

int main() {
    int steps = Parameters::STEPS;

    Mesh *mesh = new Mesh();
    mesh->InitMesh(mesh);
    std::cout << "number of cells = " << mesh->cells.size() << std::endl;

    MeshUtils::calculateEdgesNormals(*mesh);
    MeshUtils::calculateTransferVectors(*mesh);
    Initializer::initialize(*mesh);

    OutputUtils::OutputParaview(mesh, 0);

    Solver solver(mesh);
    double time = 0;
    double tau = solver.calculateTau();

    for (int i = 0; i < steps; i++) {
        std::cout << "step = " << i << " tau = " << tau << " time = " << time << std::endl;
        time += tau;

        solver.processPhase1(tau);
        solver.processPhase2(tau);
        solver.processPhase3(tau);
        solver.prepareNextStep();
        if (Parameters::HARMONIZATION_STEP > 0 && steps % Parameters::HARMONIZATION_STEP == 0) {
            solver.harmonise();
        }

        OutputUtils::OutputParaview(mesh, i + 1);
    }

    return 0;
}
