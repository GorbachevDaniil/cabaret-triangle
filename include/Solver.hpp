#ifndef Solver_hpp
#define Solver_hpp

#include "Mesh.hpp"

class Solver {
private:
    double cfl;
    Mesh *mesh;
public:
    Solver(double cfl, Mesh *mesh) : cfl(cfl), mesh(mesh) {};

    double calculateTau();
    void processPhase1(double tau);
    void processPhase2(double tau);
    void processPhase3(double tau);
    void prepareNextStep();

private:
    double calculateDivOnEdge(Edge *edge, int phase);
};

#endif
