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
    ShallowWaterSolver(double cfl, double g, Mesh *mesh) : cfl(cfl), g(g), mesh(mesh){};

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
    double calcDiv1(double h, Vector u, Vector n);
    double calcDiv2(double h, Vector u, Vector n);
    double calcDiv3(double h, Vector u, Vector n);

    double get2SqrtGH(Data *data);
    double getUx(Data *data);
    double getUy(Data *data);

    Vector calcGrad(Cell *cell, std::function<double(ShallowWaterSolver &, Data*)> getVar);

    Vector calcGradForR(long cellID, Vector n);
    Vector calcGradForQ(long cellID, Vector n);
    Vector calcGradForS(long cellID, Vector n);

    double calcInvR(double h, Vector u, double G, Vector n);
    double calcInvQ(double h, Vector u, double G, Vector n);
    double calcInvS(double h, Vector u, double G, Vector n);

    double calcLambdaR(double h, Vector u, Vector n);
    double calcLambdaQ(double h, Vector u, Vector n);
    double calcLambdaS(double h, Vector u, Vector n);

    arma::vec convertInvToInitialVariables(std::vector<arma::vec> invs);

    std::vector<arma::vec> getInvFromCellExtr(Node *node, Cell *cell, double tau);
    std::vector<arma::vec> getInvFromCellIntr(Node *node, Cell *cell, double avgH, Vector avgU, 
                                              double tau);
    
    Cell *chooseCell(
        Node *node, double avgH, Vector avgU, 
        std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda);

    Cell *chooseOppositeCell(Node *node, Cell *cell);

    double extrapolateInv(
        Cell *cell, Node *node, Vector n, double G,
        std::function<double(ShallowWaterSolver &, double, Vector, double, Vector)> calcInv,
        std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda,
        std::function<Vector(ShallowWaterSolver &, long, Vector)> calcGradForInv, 
        double tau, bool needMonotize);

    double interpolateInv(
        Cell *cell, Node *node, double avgH, Vector avgU, Vector n, double G,
        std::function<double(ShallowWaterSolver &, double, Vector, double, Vector)> calcInv,
        std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda,
        std::function<Vector(ShallowWaterSolver &, long, Vector)> calcGradForInv,
        double tau);
};

#endif
