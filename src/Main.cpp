#include "Initializer.hpp"
#include "Mesh.hpp"
#include "MeshUtils.hpp"
#include "OutputUtils.hpp"
#include "Solver.hpp"

#include <iostream>
#include <libconfig.h++>

int main() {
    libconfig::Config config;
    config.readFile("config/config.cfg");

    double cfl = config.lookup("cfl");
    int steps = config.lookup("steps");
    int writePeriod = config.lookup("write_period");
    int edgeInnerNodesNumber = config.lookup("edge_inner_nodes_number");
    bool edgeOuterNodesUsed = config.lookup("edge_outer_nodes_used");

    Mesh *mesh = new Mesh(edgeInnerNodesNumber, edgeOuterNodesUsed);
    mesh->InitMesh(mesh);
    std::cout << "number of cells = " << mesh->cells.size() << std::endl;

    MeshUtils::calculateEdgesNormals(*mesh);
    MeshUtils::calculateTransferVectors(*mesh);
    Initializer::initialize(*mesh);

    OutputUtils::OutputParaview(mesh, 0, writePeriod);

    Solver solver(cfl, mesh);
    double time = 0;
    double tau = solver.calculateTau();

    for (int i = 0; i < steps; i++) {
        std::cout << "step = " << i << " tau = " << tau << " time = " << time << std::endl;
        time += tau;

        solver.processPhase1(tau);
        solver.processPhase2(tau);
        solver.processPhase3(tau);
        solver.prepareNextStep();

        OutputUtils::OutputParaview(mesh, i + 1, writePeriod);
    }

    return 0;
}
