#include "Mesh.hpp"
#include "Solver.hpp"
#include "OutputUtils.hpp"

#include <iostream>

int main(int argc, char **argv) {
    std::cout << "Hello world, CABARET inda house" << std::endl;

    Mesh *mesh = new Mesh();
    mesh->InitMesh(mesh);

    Solver solver(mesh);

//    solver.processPhase1();
//    solver.processPhase2();
//    solver.processPhase3();

    // OutputUtils::OutputParaview(mesh);

    return 0;
}
