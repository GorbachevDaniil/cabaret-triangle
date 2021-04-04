#include "equations/shallow_water/solver.hpp"

#include <armadillo>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <limits>

double ShallowWaterSolver::calc_tau() {
    double tau = std::numeric_limits<double>::max();
    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        long centerNodeID = cell->centerNodeID;
        Data *cellData = &mesh->nodes[cell->centerNodeID].data;
        for (unsigned long edgeID : cell->edgeIDs) {
            for (unsigned long nodeID : mesh->edges[edgeID].innerNodeIDs) {
                Data *data = &mesh->nodes[nodeID].data;

                double h = 2 * (data->coords - cellData->coords).length();
                double lambdaR = calcLambdaR(mesh->s0[centerNodeID][0], mesh->v0[centerNodeID][0],
                                             cell->nodeToTransferVector[nodeID]);
                double nodeTauR = h / std::abs(lambdaR);
                if (tau > nodeTauR) {
                    tau = nodeTauR;
                }
                double lambdaQ = calcLambdaQ(mesh->s0[centerNodeID][0], mesh->v0[centerNodeID][0],
                                             cell->nodeToTransferVector[nodeID]);
                double nodeTauQ = h / std::abs(lambdaQ);
                if (tau > nodeTauQ) {
                    tau = nodeTauQ;
                }
            }
        }
    }
    assert(tau != std::numeric_limits<double>::max());
    return cfl * tau;
}

double ShallowWaterSolver::calcDiv1(double h, Vector u, Vector n) {
    return (u * h) * n;
}

double ShallowWaterSolver::calcDiv2(double h, Vector u, Vector n) {
    Vector F = Vector(pow(u.x, 2) + g * h / 2, u.x * u.y);
    return (F * h) * n;
}

double ShallowWaterSolver::calcDiv3(double h, Vector u, Vector n) {
    Vector F = Vector(u.x * u.y, pow(u.y, 2) + g * h / 2);
    return (F * h) * n;
}

double ShallowWaterSolver::calcIntegral(std::vector<double> values, double length) {
    if (values.size() == 4) {
        return (values[0] + 3 * values[1] + 3 * values[2] + values[3]) / 8 * length;
    }
    if (values.size() == 2) {
        return (values[0] + values[1]) / 2 * length;
    }
    assert(false && "such integration not supported");
}

void ShallowWaterSolver::process_phase_1(double tau) {
    for (unsigned long i = 0; i < mesh->edges.size(); i++) {
        Edge *edge = &mesh->edges[i];

        std::vector<double> div1Edge;
        std::vector<double> div2Edge;
        std::vector<double> div3Edge;

        Vector n = edge->normal;

        for (unsigned long usedNodeID : edge->usedNodeIDs) {
            double h = mesh->s0[usedNodeID][0];
            Vector u = mesh->v0[usedNodeID][0];

            div1Edge.push_back(calcDiv1(h, u, n));
            div2Edge.push_back(calcDiv2(h, u, n));
            div3Edge.push_back(calcDiv3(h, u, n));
        }

        mesh->edgeIDToDivs0[edge->ID][0] = calcIntegral(div1Edge, edge->length);
        mesh->edgeIDToDivs0[edge->ID][1] = calcIntegral(div2Edge, edge->length);
        mesh->edgeIDToDivs0[edge->ID][2] = calcIntegral(div3Edge, edge->length);
    }

    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];

        double div1 = 0;
        double div2 = 0;
        double div3 = 0;
        for (long edgeID : cell->edgeIDs) {
            int coef = cell->edgeToNormalDir[edgeID];

            div1 += coef * mesh->edgeIDToDivs0[edgeID][0];
            div2 += coef * mesh->edgeIDToDivs0[edgeID][1];
            div3 += coef * mesh->edgeIDToDivs0[edgeID][2];
        }
        div1 /= cell->volume;
        div2 /= cell->volume;
        div3 /= cell->volume;

        long centerNodeID = cell->centerNodeID;

        double h = mesh->s0[centerNodeID][0];
        Vector u = mesh->v0[centerNodeID][0];
        double newH = h - tau * div1 / 2;
        double newUX = (h * u.x - tau * div2 / 2) / newH;
        double newUY = (h * u.y - tau * div3 / 2) / newH;

        mesh->s1[centerNodeID][0] = newH;
        mesh->v1[centerNodeID][0].x = newUX;
        mesh->v1[centerNodeID][0].y = newUY;
    }
}

double ShallowWaterSolver::calcInvR(double G, double h, Vector u, Vector n) {
    return u * n + G * h;
}

double ShallowWaterSolver::calcInvQ(double G, double h, Vector u, Vector n) {
    return u * n - G * h;
}

double ShallowWaterSolver::calcInvS(double /*G*/, double /*h*/, Vector u, Vector n) {
    return -u.x * n.y + u.y * n.x;
}

double ShallowWaterSolver::calcLambdaR(double h, Vector u, Vector n) {
    return u * n + sqrt(g * h);
}

double ShallowWaterSolver::calcLambdaQ(double h, Vector u, Vector n) {
    return u * n - sqrt(g * h);
}

double ShallowWaterSolver::calcLambdaS(double /*h*/, Vector u, Vector n) {
    return u * n;
}

double ShallowWaterSolver::calcQR(Cell *cell, double dirDerivH,
                                  double dirDerivUx, double dirDerivUy,
                                  Vector n) {
    Vector m = Vector(-n.y, n.x);
    long center = cell->centerNodeID;
    double Q = 0;
    double centerDirVelocity = m * mesh->v1[center][0];
    double centerH = mesh->s1[center][0];
    Q += sqrt(g / centerH) * centerDirVelocity * dirDerivH;
    Q += (m.x * sqrt(g * centerH) + n.x * centerDirVelocity) * dirDerivUx;
    Q += (m.y * sqrt(g * centerH) + n.y * centerDirVelocity) * dirDerivUy;
    return -Q;
}

double ShallowWaterSolver::calcQQ(Cell *cell, double dirDerivH,
                                  double dirDerivUx, double dirDerivUy,
                                  Vector n) {
    Vector m = Vector(-n.y, n.x);
    long center = cell->centerNodeID;
    double Q = 0;
    double centerDirVelocity = m * mesh->v1[center][0];
    double centerH = mesh->s1[center][0];
    Q += -sqrt(g / centerH) * centerDirVelocity * dirDerivH;
    Q += (-m.x * sqrt(g * centerH) + n.x * centerDirVelocity) * dirDerivUx;
    Q += (-m.y * sqrt(g * centerH) + n.y * centerDirVelocity) * dirDerivUy;
    return -Q;
}

double ShallowWaterSolver::calcQS(Cell *cell, double dirDerivH,
                                  double dirDerivUx, double dirDerivUy,
                                  Vector n) {
    Vector m = Vector(-n.y, n.x);
    long center = cell->centerNodeID;
    double Q = 0;
    Q += g * dirDerivH;
    Q += (m * mesh->v1[center][0]) * (m.x * dirDerivUx + m.y * dirDerivUy);
    return -Q;
}

double monotize(double inv2, double inv0, double invCenter0, double invOpposite0, double invCenter1,
                double lambda, double tau, double /*h*/, double dirDerivative) {

    double Q = (invCenter1 - invCenter0) / (tau / 2) + lambda * dirDerivative;
    // invCenter0 = (6 * invCenter0 / h - inv0 - invOpposite0) / 4;
    double min = std::min(std::min(inv0, invOpposite0), invCenter0);
    double max = std::max(std::max(inv0, invOpposite0), invCenter0);
    min += tau * Q;
    max += tau * Q;

    inv2 = std::max(inv2, min);
    inv2 = std::min(inv2, max);

    return inv2;
}

double monotize(double inv2, double inv0, double invCenter0, double invOpposite0,
                double /*invCenter1*/, double tau, double Q) {
    double min = std::min(std::min(inv0, invOpposite0), invCenter0);
    double max = std::max(std::max(inv0, invOpposite0), invCenter0);

    min += tau * Q;
    max += tau * Q;

    inv2 = std::max(inv2, min);
    inv2 = std::min(inv2, max);

    return inv2;
}

double ShallowWaterSolver::extrapolateInv(
    Cell* cell, Edge* /*edge*/, Node* node, Vector n, double G,
    std::function<double(ShallowWaterSolver &, double, double, Vector, Vector)> calcInv,
    std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda,
    std::function<double(ShallowWaterSolver &, Cell *, double, double, double, Vector)> /*calcQ*/,
    double tau, bool needMonotize) {
    assert(cell->nodeIDToOppositeNodeID[node->ID] != -1);

    long considered = node->ID;
    long center = cell->centerNodeID;
    long opposite = cell->nodeIDToOppositeNodeID[node->ID];
    double h = (mesh->nodes[considered].data.coords - mesh->nodes[opposite].data.coords).length();

    double invCenter1 = calcInv(*this, G, mesh->s1[center][0], mesh->v1[center][0], n);
    double invOpposite0 = calcInv(*this, G, mesh->s0[opposite][0], mesh->v0[opposite][0], n);
    double coefPan = 0.0;
    double inv2 = (2 * invCenter1 - (1 - coefPan) * invOpposite0) / (1 + coefPan);
    if (needMonotize) {
        double inv0 = calcInv(*this, G, mesh->s0[considered][0], mesh->v0[considered][0], n);
        double invCenter0 = calcInv(*this, G, mesh->s0[center][0], mesh->v0[center][0], n);
        double lambda = calcLambda(*this, mesh->s1[center][0], mesh->v1[center][0], n);

        Vector transfer = cell->nodeToTransferVector[considered];
        Vector grad = Vector(0, 0);
        for (long edgeID : cell->edgeIDs) {
            Edge *tmpEdge = &mesh->edges[edgeID];
            std::vector<double> variables;
            for (unsigned long usedNodeID : tmpEdge->usedNodeIDs) {
                double inv = calcInv(*this, G,
                                     mesh->s0[usedNodeID][0],
                                     mesh->v0[usedNodeID][0],
                                     n);
                variables.push_back(inv);
            }
            Vector edgeN = tmpEdge->normal * cell->edgeToNormalDir[edgeID];
            grad = grad + edgeN * calcIntegral(variables, tmpEdge->length);
        }
        grad = grad / cell->volume;

        inv2 = monotize(inv2, inv0, invCenter0, invOpposite0, invCenter1,
                        lambda, tau, h, transfer * grad);


        // Vector transfer = cell->nodeToTransferVector[considered];
        // double coef = cell->edgeToNormalDir[edge->ID];

        // double h = (mesh->nodes[considered].data.coords - mesh->nodes[opposite].data.coords).length();
        // double derivative = coef * (inv0 - invOpposite0) / h;
        // inv2 = monotize(inv2, inv0, invCenter0, invOpposite0, invCenter1,
        //                 lambda, tau, derivative);

        // Vector transfer = cell->nodeToTransferVector[considered];
        // Vector oppositeTransfer = Vector(-transfer.y, transfer.x);
        // Vector gradH = Vector(0, 0);
        // Vector gradUx = Vector(0, 0);
        // Vector gradUy = Vector(0, 0);
        // for (long edgeID : cell->edgeIDs) {
        //     Edge *tmpEdge = &mesh->edges[edgeID];
        //     std::vector<double> h;
        //     std::vector<double> ux;
        //     std::vector<double> uy;
        //     for (unsigned long usedNodeID : tmpEdge->usedNodeIDs) {
        //         h.push_back(mesh->s0[usedNodeID][0]);
        //         ux.push_back(mesh->v0[usedNodeID][0].x);
        //         uy.push_back(mesh->v0[usedNodeID][0].y);
        //     }
        //     Vector edgeN = tmpEdge->normal * cell->edgeToNormalDir[edgeID];
        //     gradH = gradH + edgeN * calcIntegral(h, tmpEdge->length);
        //     gradUx = gradUx + edgeN * calcIntegral(ux, tmpEdge->length);
        //     gradUy = gradUy + edgeN * calcIntegral(uy, tmpEdge->length);
        // }
        // gradH = gradH / cell->volume;
        // gradUx = gradUx / cell->volume;
        // gradUy = gradUy / cell->volume;
        // double dirDerivativeH = gradH * oppositeTransfer;
        // double dirDerivativeUx = gradUx * oppositeTransfer;
        // double dirDerivativeUy = gradUy * oppositeTransfer;

        // double Q = calcQ(*this, cell, dirDerivativeH, dirDerivativeUx, dirDerivativeUy, transfer);
        // inv2 = monotize(inv2, inv0, invCenter0, invOpposite0, invCenter1,
        //                 tau, Q);
    }

    return inv2;
}

arma::vec ShallowWaterSolver::convertInvToInitialVariables(std::vector<arma::vec> invs) {
    assert(invs.size() == 3);

    arma::mat A = arma::mat(3, 3);
    arma::vec b = arma::vec(3);

    A(0, 0) = invs[0][0];
    A(0, 1) = invs[0][1];
    A(0, 2) = invs[0][2];
    b[0] = invs[0][3];

    A(1, 0) = invs[1][0];
    A(1, 1) = invs[1][1];
    A(1, 2) = invs[1][2];
    b[1] = invs[1][3];

    A(2, 0) = invs[2][0];
    A(2, 1) = invs[2][1];
    A(2, 2) = invs[2][2];
    b[2] = invs[2][3];

    if (arma::rank(A) != 3) {
        std::cout << A << b << std::endl;
        assert(false);
    }
    arma::vec initialVariables = arma::solve(A, b);
    return initialVariables;
}

std::vector<arma::vec> ShallowWaterSolver::getInvFromCellExtr(Node *node, Edge *edge, Cell *cell,
                                                              double tau) {
    bool monotize = true;

    std::vector<arma::vec> invs;

    long centerNodeID = cell->centerNodeID;
    double h = mesh->s1[centerNodeID][0];
    Vector u = mesh->v1[centerNodeID][0];
    double G = sqrt(g / h);

    double coef = 1;
    Vector n = edge->normal * cell->edgeToNormalDir[edge->ID];

    // double coef = cell->edgeToNormalDir[edge->ID];
    // Vector n = edge->normal;

    double lambdaR = coef * calcLambdaR(h, u, n);
    if (lambdaR > 0) {
        double R = extrapolateInv(cell, edge, node, n, G,
                                  &ShallowWaterSolver::calcInvR,
                                  &ShallowWaterSolver::calcLambdaR,
                                  &ShallowWaterSolver::calcQR,
                                  tau, monotize);
        arma::vec inv(4);
        inv[0] = n.x;
        inv[1] = n.y;
        inv[2] = G;
        inv[3] = R;
        invs.push_back(inv);
    }

    double lambdaQ = coef * calcLambdaQ(h, u, n);
    if (lambdaQ > 0) {
        double Q = extrapolateInv(cell, edge, node, n, G,
                                  &ShallowWaterSolver::calcInvQ,
                                  &ShallowWaterSolver::calcLambdaQ,
                                  &ShallowWaterSolver::calcQQ,
                                  tau, monotize);
        arma::vec inv(4);
        inv[0] = n.x;
        inv[1] = n.y;
        inv[2] = -G;
        inv[3] = Q;
        invs.push_back(inv);
    }

    double lambdaS = coef * calcLambdaS(h, u, n);
    double S = extrapolateInv(cell, edge, node, n, G,
                              &ShallowWaterSolver::calcInvS,
                              &ShallowWaterSolver::calcLambdaS,
                              &ShallowWaterSolver::calcQS,
                              tau, monotize);
    arma::vec inv(5);
    inv[0] = -n.y;
    inv[1] = n.x;
    inv[2] = 0;
    inv[3] = S;
    inv[4] = lambdaS;
    invs.push_back(inv);

    assert(invs.size() == 2);
    return invs;
}

void ShallowWaterSolver::processPhase2BoundEdge(Edge *edge, double tau) {
    Cell *cell = &mesh->cells[edge->cellIDs[0]];
    long centerNodeID = cell->centerNodeID;

    for (unsigned long nodeID : edge->innerNodeIDs) {
        Node *node = &mesh->nodes[nodeID];

        std::vector<arma::vec> invs = getInvFromCellExtr(node, edge, cell, tau);

        // --- no-slip condition ---
        // mesh->v2[nodeID][0].x = 0;
        // mesh->v2[nodeID][0].y = 0;
        // mesh->s2[nodeID][0] = mesh->s1[centerNodeID][0];
        // --- ---

        // --- slip condition ---
        Vector edgeN = edge->normal * cell->edgeToNormalDir[edge->ID];
        arma::vec inv(4);
        inv[0] = edgeN.x;
        inv[1] = edgeN.y;
        inv[2] = 0;
        inv[3] = 0;
        invs.push_back(inv);

        if (invs[1][4] <= 0) {
            invs[1][3] = calcInvS(0, 0, mesh->v1[centerNodeID][0], edgeN);
        }

        arma::vec initialValues = convertInvToInitialVariables(invs);

        mesh->v2[nodeID][0].x = initialValues[0];
        mesh->v2[nodeID][0].y = initialValues[1];
        mesh->s2[nodeID][0] = initialValues[2];
        // --- ---

        // --- flow condition ---
        // Vector edgeN = edge->normal * cell->edgeToNormalDir[edge->ID];
        // arma::vec inv(4);
        // inv[0] = -edgeN.x;
        // inv[1] = -edgeN.y;
        // inv[2] = 2;
        // inv[3] = 2 * sqrt(g * 1);
        // invs.push_back(inv);

        // if (invs[1][4] <= 0) {
        //     Vector n = cell->nodeToTransferVector[nodeID];
        //     invs[1][3] = calcInvS(0, 0, mesh->v1[centerNodeID][0], n);
        // }

        // arma::vec initialValues = convertInvToInitialVariables(invs);

        // mesh->v2[nodeID][0].x = initialValues[0];
        // mesh->v2[nodeID][0].y = initialValues[1];
        // mesh->s2[nodeID][0] = initialValues[2];
        // --- ---
    }
}

void ShallowWaterSolver::processPhase2InnerEdge(Edge *edge, double tau) {
    Cell *cell1 = &mesh->cells[edge->cellIDs[0]];
    Cell *cell2 = &mesh->cells[edge->cellIDs[1]];

    for (unsigned long nodeID : edge->innerNodeIDs) {
        Node *node = &mesh->nodes[nodeID];

        std::vector<arma::vec> invs1 = getInvFromCellExtr(node, edge, cell1, tau);
        std::vector<arma::vec> invs2 = getInvFromCellExtr(node, edge, cell2, tau);
        std::vector<arma::vec> invs;

        invs.push_back(invs1[0]);
        invs.push_back(invs2[0]);
        if (invs1[1][4] > invs2[1][4]) {
            invs.push_back(invs1[1]);
        } else {
            invs.push_back(invs2[1]);
        }
        arma::vec initialValues = convertInvToInitialVariables(invs);

        mesh->v2[nodeID][0].x = initialValues[0];
        mesh->v2[nodeID][0].y = initialValues[1];
        mesh->s2[nodeID][0] = initialValues[2];
    }
}

void ShallowWaterSolver::process_phase_2(double tau) {
    for (unsigned long i = 0; i < mesh->edges.size(); i++) {
        Edge *edge = &mesh->edges[i];
        if (edge->boundEdge) {
            processPhase2BoundEdge(edge, tau);
        } else {
            processPhase2InnerEdge(edge, tau);
        }
    }
}

void ShallowWaterSolver::process_phase_3(double tau) {
    for (unsigned long i = 0; i < mesh->edges.size(); i++) {
        Edge *edge = &mesh->edges[i];

        std::vector<double> div1Edge;
        std::vector<double> div2Edge;
        std::vector<double> div3Edge;

        Vector n = edge->normal;

        for (unsigned long usedNodeID : edge->usedNodeIDs) {
            double h = mesh->s2[usedNodeID][0];
            Vector u = mesh->v2[usedNodeID][0];

            div1Edge.push_back(calcDiv1(h, u, n));
            div2Edge.push_back(calcDiv2(h, u, n));
            div3Edge.push_back(calcDiv3(h, u, n));
        }

        mesh->edgeIDToDivs2[edge->ID][0] = calcIntegral(div1Edge, edge->length);
        mesh->edgeIDToDivs2[edge->ID][1] = calcIntegral(div2Edge, edge->length);
        mesh->edgeIDToDivs2[edge->ID][2] = calcIntegral(div3Edge, edge->length);
    }

    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];

        double div1 = 0;
        double div2 = 0;
        double div3 = 0;
        for (long edgeID : cell->edgeIDs) {
            int coef = cell->edgeToNormalDir[edgeID];

            div1 += coef * (mesh->edgeIDToDivs2[edgeID][0] + mesh->edgeIDToDivs0[edgeID][0]);
            div2 += coef * (mesh->edgeIDToDivs2[edgeID][1] + mesh->edgeIDToDivs0[edgeID][1]);
            div3 += coef * (mesh->edgeIDToDivs2[edgeID][2] + mesh->edgeIDToDivs0[edgeID][2]);
        }
        div1 /= cell->volume;
        div2 /= cell->volume;
        div3 /= cell->volume;

        long centerNodeID = cell->centerNodeID;

        double h = mesh->s0[centerNodeID][0];
        Vector u = mesh->v0[centerNodeID][0];
        double newH = h - tau * div1 / 2;
        double newUX = (h * u.x - tau * div2 / 2) / newH;
        double newUY = (h * u.y - tau * div3 / 2) / newH;

        mesh->s2[centerNodeID][0] = newH;
        mesh->v2[centerNodeID][0].x = newUX;
        mesh->v2[centerNodeID][0].y = newUY;
    }
}

void ShallowWaterSolver::prepare_next_step() {
    for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
        for (unsigned long j = 0; j < mesh->s0[i].size(); j++) {
            mesh->s0[i][j] = mesh->s2[i][j];
        }
        for (unsigned long j = 0; j < mesh->v0[i].size(); j++) {
            mesh->v0[i][j] = mesh->v2[i][j];
        }
    }
}