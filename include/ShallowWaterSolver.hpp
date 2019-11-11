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

    void processPhase2BoundNode(Node *node, double tau);
    void processPhase2InnerNode(Node *node, double tau);

   private:
    double calcIntegral(std::vector<double> values, double length);

    double calcDiv1(double h, Vector u, Vector n);
    double calcDiv2(double h, Vector u, Vector n);
    double calcDiv3(double h, Vector u, Vector n);

    double get2SqrtGH(long nodeID);
    double getUx(long nodeID);
    double getUy(long nodeID);

    Vector calcGrad(Cell *cell, std::function<double(ShallowWaterSolver &, long)> getVar);

    Vector calcGradForR(long cellID, Vector n);
    Vector calcGradForQ(long cellID, Vector n);
    Vector calcGradForS(long cellID, Vector n);

    double calcInvR(double h, Vector u, Vector n);
    double calcInvQ(double h, Vector u, Vector n);
    double calcInvS(double h, Vector u, Vector n);

    double calcLambdaR(double h, Vector u, Vector n);
    double calcLambdaQ(double h, Vector u, Vector n);
    double calcLambdaS(double h, Vector u, Vector n);

    arma::vec convertInvToInitialVariables(std::vector<arma::vec> invs);

    std::vector<arma::vec> getInvFromCellExtr(Node *node, Edge *edge, Cell *cell, double tau);
    std::vector<arma::vec> getInvFromCellIntr(Node *node, Cell *cell, double avgH, Vector avgU,
                                              double tau);
    std::vector<arma::vec> getInvFromEdge(Node *node, Edge *edge, double avgH, Vector avgU,
                                          double tau);

    double extrapolateInv(
        Cell *cell, Node *node, Vector n,
        std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcInv,
        std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda,
        std::function<Vector(ShallowWaterSolver &, long, Vector)> calcGradForInv, double tau,
        bool needMonotize);

    // --- calculate new values by interpolation in cell ---
    Cell *chooseCell(
        Node *node, double avgH, Vector avgU,
        std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda);
    Cell *chooseOppositeCell(Node *node, Cell *cell);
    double interpolateInv(
        Cell *cell, Node *node, double avgH, Vector avgU, Vector n,
        std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcInv,
        std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda,
        std::function<Vector(ShallowWaterSolver &, long, Vector)> calcGradForInv, double tau);
    // -------

    // --- calculate new values by using values from edges ---
    Edge *chooseEdge(
        Node *node, double avgH, Vector avgU,
        std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda);
    Edge *chooseOppositeEdge(Node *node, Edge *edge);
    double extrapolateInv(
        Edge *edge, Node *node, Vector n,
        std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcInv,
        std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda,
        std::function<Vector(ShallowWaterSolver &, long, Vector)> calcGradForInv, double tau,
        bool needMonotize);
    // -------
};

#endif
