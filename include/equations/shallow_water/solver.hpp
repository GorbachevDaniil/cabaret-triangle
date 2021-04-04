#ifndef ShallowWaterSolver_hpp
#define ShallowWaterSolver_hpp

#include <functional>

#include "../virtual_solver.hpp"

class ShallowWaterSolver : public Solver {
public:
    ShallowWaterSolver(double cfl, double g, Mesh *mesh) {
        this->cfl = cfl;
        this->g = g;
        this->mesh = mesh;
    };

    double calc_tau();
    void process_phase_1(double tau);
    void process_phase_2(double tau);
    void process_phase_3(double tau);
    void prepare_next_step();

private:
    void processPhase2BoundEdge(Edge *edge, double tau);
    void processPhase2InnerEdge(Edge *edge, double tau);

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
        std::function<double(ShallowWaterSolver &, Cell *, double, double, double, Vector)> calcQ,
        double tau, bool needMonotize);

private:
    double cfl;
    double g;

    Mesh *mesh;
};

#endif
