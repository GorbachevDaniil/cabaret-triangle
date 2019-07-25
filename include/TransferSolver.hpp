#ifndef TransferSolver_hpp
#define TransferSolver_hpp

#include "AbstractSolver.hpp"

class TransferSolver : public AbstractSolver {
   private:
    double cfl;
    Mesh *mesh;

   public:
    TransferSolver(double cfl, Mesh *mesh) {
        this->cfl = cfl;
        this->mesh = mesh;
    };

    double calcTau();
    void processPhase1(double tau);
    void processPhase2(double tau);
    void processPhase3(double tau);
    void prepareNextStep();

    void processPhase2BoundEdge(Edge *edge, double tau);
    void processPhase2InnerEdge(Edge *edge, double tau);

   private:
    double calcDivOnEdge(Edge *edge, int phase);
    double getNewInvariantValue(Data *data, Data *centerData, Data *oppositeData, Vector transfer,
                                double tau, bool monotize);
};

#endif
