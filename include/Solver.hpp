#ifndef Solver_hpp
#define Solver_hpp

#include "Mesh.hpp"

class Solver {
private:
    Mesh *mesh;
public:
    Solver(Mesh *mesh) : mesh(mesh) {};

    double calculateTau();
    void processPhase1(double tau);
    void processPhase2(double tau);
    void processPhase3(double tau);
    void prepareNextStep();
    void harmonise();

private:
    double calculateDivOnEdge(Edge *edge, int phase);
};

#endif
