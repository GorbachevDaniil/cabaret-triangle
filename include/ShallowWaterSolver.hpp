#ifndef ShallowWaterSolver_hpp
#define ShallowWaterSolver_hpp

#include "AbstractSolver.hpp"

#include <functional>

class ShallowWaterSolver : public AbstractSolver {
   private:
    double cfl;
    double g;

    Mesh *mesh;

   public:
    ShallowWaterSolver(double cfl, double g, Mesh *mesh) {
        this->cfl = cfl;
        this->g = g;
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
    double calcIntegral(std::vector<double> values, double length);

    double calcDiv1(double h, Vector u, Vector n);
    double calcDiv2(double h, Vector u, Vector n);
    double calcDiv3(double h, Vector u, Vector n);

    double calcInvR(double G, double h, Vector u, Vector n);
    double calcInvQ(double G, double h, Vector u, Vector n);
    double calcInvS(double G, double h, Vector u, Vector n);

    double calcLambdaR(double h, Vector u, Vector n);
    double calcLambdaQ(double h, Vector u, Vector n);
    double calcLambdaS(double h, Vector u, Vector n);

    double calcQR(Cell *cell, double dirDerivH, double dirDerivUx, double dirDerivUy, Vector n);
    double calcQQ(Cell *cell, double dirDerivH, double dirDerivUx, double dirDerivUy, Vector n);
    double calcQS(Cell *cell, double dirDerivH, double dirDerivUx, double dirDerivUy, Vector n);

    arma::vec convertInvToInitialVariables(std::vector<arma::vec> invs);

    std::vector<arma::vec> getInvFromCellExtr(Node *node, Edge *edge, Cell *cell, double tau);

    double extrapolateInv(
        Cell *cell, Edge *edge, Node *node, Vector n, double G,
        std::function<double(ShallowWaterSolver &, double, double, Vector, Vector)> calcInv,
        std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda,
        std::function<double(ShallowWaterSolver &, Cell*, double, double, double, Vector)> calcQ,
        double tau, bool needMonotize);
};

#endif
