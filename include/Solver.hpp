#ifndef Solver_hpp
#define Solver_hpp

#include "Mesh.hpp"

class Solver {
private:
    Mesh *mesh;
public:
    Solver(Mesh *mesh) : mesh(mesh) {};

    void processPhase1();
    void processPhase2();
    void processPhase3();
    void prepareNextStep();
};

#endif
