#include "ShallowWaterSolver.hpp"

#include <armadillo>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <limits>

double ShallowWaterSolver::calcTau() {
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

double ShallowWaterSolver::calcDiv1(double h, Vector u, Vector n) { return u * n * h; }

double ShallowWaterSolver::calcDiv2(double h, Vector u, Vector n) {
    Vector F = Vector(pow(u.x, 2) + g * h / 2, u.x * u.y);
    return F * n * h;
}

double ShallowWaterSolver::calcDiv3(double h, Vector u, Vector n) {
    Vector F = Vector(u.x * u.y, pow(u.y, 2) + g * h / 2);
    return F * n * h;
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

void ShallowWaterSolver::processPhase1(double tau) {
    for (unsigned long i = 0; i < mesh->edges.size(); i++) {
        Edge *edge = &mesh->edges[i];

        std::vector<double> div1Edge;
        std::vector<double> div2Edge;
        std::vector<double> div3Edge;

        Vector n = edge->normal;

        for (unsigned long usedNodeID : edge->getUsedNodes(*mesh)) {
            double h = mesh->s0[usedNodeID][0];
            Vector u = mesh->v0[usedNodeID][0];

            div1Edge.push_back(calcDiv1(h, u, n));
            div2Edge.push_back(calcDiv2(h, u, n));
            div3Edge.push_back(calcDiv3(h, u, n));
        }

        mesh->edgeIDToDivs[edge->ID][0] = calcIntegral(div1Edge, edge->length);
        mesh->edgeIDToDivs[edge->ID][1] = calcIntegral(div2Edge, edge->length);
        mesh->edgeIDToDivs[edge->ID][2] = calcIntegral(div3Edge, edge->length);
    }

    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];

        double div1 = 0;
        double div2 = 0;
        double div3 = 0;
        for (long edgeID : cell->edgeIDs) {
            int coef = cell->edgeToNormalDir[edgeID];

            div1 += coef * mesh->edgeIDToDivs[edgeID][0];
            div2 += coef * mesh->edgeIDToDivs[edgeID][1];
            div3 += coef * mesh->edgeIDToDivs[edgeID][2];
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

double ShallowWaterSolver::calcInvR(double h, Vector u, Vector n) {
    return u * n + 2 * sqrt(g * h);
}

double ShallowWaterSolver::calcInvQ(double h, Vector u, Vector n) {
    return u * n - 2 * sqrt(g * h);
}

double ShallowWaterSolver::calcInvS(double h, Vector u, Vector n) { return -u.x * n.y + u.y * n.x; }

double ShallowWaterSolver::calcLambdaR(double h, Vector u, Vector n) { return u * n + sqrt(g * h); }

double ShallowWaterSolver::calcLambdaQ(double h, Vector u, Vector n) { return u * n - sqrt(g * h); }

double ShallowWaterSolver::calcLambdaS(double h, Vector u, Vector n) { return u * n; }

double monotize(double inv2, double inv0, double invCenter0, double invOpposite0, double invCenter1,
                double lambda, double tau, double dirDerivative, bool withQ) {
    double min = std::min(std::min(inv0, invCenter0), invOpposite0);
    double max = std::max(std::max(inv0, invCenter0), invOpposite0);

    if (withQ) {
        double Q = (invCenter1 - invCenter0) / (tau / 2) + lambda * dirDerivative;
        min += tau * Q;
        max += tau * Q;
    }

    inv2 = std::max(inv2, min);
    inv2 = std::min(inv2, max);

    return inv2;
}

double ShallowWaterSolver::get2SqrtGH(long nodeID) { return 2 * sqrt(g * mesh->s0[nodeID][0]); }

double ShallowWaterSolver::getUx(long nodeID) { return mesh->v0[nodeID][0].x; }

double ShallowWaterSolver::getUy(long nodeID) { return mesh->v0[nodeID][0].y; }

Vector ShallowWaterSolver::calcGrad(Cell *cell,
                                    std::function<double(ShallowWaterSolver &, long)> getVar) {
    Vector grad = Vector(0, 0);

    for (long edgeID : cell->edgeIDs) {
        Edge *edge = &mesh->edges[edgeID];
        std::vector<double> variables;
        for (unsigned long usedNodeID : edge->getUsedNodes(*mesh)) {
            variables.push_back(getVar(*this, usedNodeID));
        }
        Vector edgeN = edge->normal * cell->edgeToNormalDir[edgeID];
        grad = grad + edgeN * calcIntegral(variables, edge->length);
    }
    grad = grad / cell->volume;

    return grad;
}

Vector ShallowWaterSolver::calcGradForR(long cellID, Vector n) {
    Vector grad = mesh->cellIDToGrads[cellID][0];

    grad.x += n.x * mesh->cellIDToGrads[cellID][1].x + n.y * mesh->cellIDToGrads[cellID][2].x;
    grad.y += n.x * mesh->cellIDToGrads[cellID][1].y + n.y * mesh->cellIDToGrads[cellID][2].y;

    return grad;
}

Vector ShallowWaterSolver::calcGradForQ(long cellID, Vector n) {
    Vector grad = mesh->cellIDToGrads[cellID][0] * -1;

    grad.x += n.x * mesh->cellIDToGrads[cellID][1].x + n.y * mesh->cellIDToGrads[cellID][2].x;
    grad.y += n.x * mesh->cellIDToGrads[cellID][1].y + n.y * mesh->cellIDToGrads[cellID][2].y;

    return grad;
}

Vector ShallowWaterSolver::calcGradForS(long cellID, Vector n) {
    Vector grad = Vector(0, 0);

    grad.x += -n.y * mesh->cellIDToGrads[cellID][1].x + n.x * mesh->cellIDToGrads[cellID][2].x;
    grad.y += -n.y * mesh->cellIDToGrads[cellID][1].y + n.x * mesh->cellIDToGrads[cellID][2].y;

    return grad;
}

double ShallowWaterSolver::extrapolateInv(
    Cell *cell, Node *node, Vector n,
    std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcInv,
    std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda,
    std::function<Vector(ShallowWaterSolver &, long, Vector)> calcGradForInv, double tau,
    bool needMonotize) {
    assert(cell->nodeIDToOppositeNodeID[node->ID] != -1);

    long nodeID = node->ID;
    long centerNodeID = cell->centerNodeID;
    long oppositeNodeID = cell->nodeIDToOppositeNodeID[node->ID];

    double invCenter1 = calcInv(*this, mesh->s1[centerNodeID][0], mesh->v1[centerNodeID][0], n);
    double invOpposite0 = calcInv(*this, mesh->s0[oppositeNodeID][0], mesh->v0[oppositeNodeID][0], n);
    double coefPan = 0.0;
    double inv2 = (2 * invCenter1 - (1 - coefPan) * invOpposite0) / (1 + coefPan);
    if (needMonotize) {
        double inv0 = calcInv(*this, mesh->s0[nodeID][0], mesh->v0[nodeID][0], n);
        double invCenter0 = calcInv(*this, mesh->s0[centerNodeID][0], mesh->v0[centerNodeID][0], n);

        double lambda = calcLambda(*this, mesh->s1[centerNodeID][0], mesh->v1[centerNodeID][0], n);
        Vector transfer = cell->nodeToTransferVector[nodeID];
        Vector grad = calcGradForInv(*this, cell->ID, transfer);
        inv2 = monotize(inv2, inv0, invCenter0, invOpposite0, invCenter1, lambda, tau,
                        transfer * grad, true);
    }

    return inv2;
}

double interpolate(double x, double y, arma::vec a) {
    double value = 0;
    value += a[0];
    value += a[1] * x;
    value += a[2] * y;
    value += a[3] * x * x;
    value += a[4] * y * y;
    value += a[5] * x * y;
    value += a[6] * x * x * x;
    value += a[7] * y * y * y;
    value += a[8] * x * y * x;
    value += a[9] * x * y * y;
    return value;
}

double ShallowWaterSolver::interpolateInv(
    Cell *cell, Node *node, double avgH, Vector avgU, Vector n,
    std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcInv,
    std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda,
    std::function<Vector(ShallowWaterSolver &, long, Vector)> calcGradForInv, double tau) {
    arma::vec invs(10);
    int pos = 0;
    for (long nodeID : cell->nodeIDs) {
        invs[pos] = calcInv(*this, mesh->s0[nodeID][0], mesh->v0[nodeID][0], n);
        pos++;
    }
    for (long edgeID : cell->edgeIDs) {
        for (long nodeID : mesh->edges[edgeID].innerNodeIDs) {
            invs[pos] = calcInv(*this, mesh->s0[nodeID][0], mesh->v0[nodeID][0], n);
            pos++;
        }
    }
    long centerNodeID = cell->centerNodeID;
    invs[pos] = calcInv(*this, mesh->s0[centerNodeID][0], mesh->v0[centerNodeID][0], n);
    arma::vec inv = cell->interpolationMat * invs;

    Data *data = &mesh->nodes[node->ID].data;
    double coef = -calcLambda(*this, avgH, avgU, n) * tau;
    Vector intersectCoords = n * coef + data->coords;
    double intersectX = intersectCoords.x;
    double intersectY = intersectCoords.y;

    double inv2 = interpolate(intersectX, intersectY, inv);

    double invCenter0 = calcInv(*this, mesh->s0[centerNodeID][0], mesh->v0[centerNodeID][0], n);
    double invCenter1 = calcInv(*this, mesh->s1[centerNodeID][0], mesh->v1[centerNodeID][0], n);
    double lambda = calcLambda(*this, mesh->s1[centerNodeID][0], mesh->v1[centerNodeID][0], n);

    Vector grad = calcGradForInv(*this, cell->ID, n);
    inv2 += tau * ((invCenter1 - invCenter0) / (tau / 2) + lambda * (grad * n));

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
    initialVariables[2] = pow(initialVariables[2], 2) / g;
    return initialVariables;
}

std::vector<arma::vec> ShallowWaterSolver::getInvFromCellExtr(Node *node, Edge *edge, Cell *cell,
                                                              double tau) {
    std::vector<arma::vec> invs;

    long centerNodeID = cell->centerNodeID;
    double h = mesh->s1[centerNodeID][0];
    Vector u = mesh->v1[centerNodeID][0];

    // Vector n = cell->nodeToTransferVector[node->ID];
    Vector n = edge->normal * cell->edgeToNormalDir[edge->ID];

    double lambdaR = calcLambdaR(h, u, n);
    if (lambdaR > 0) {
        double R = extrapolateInv(cell, node, n,
                                  &ShallowWaterSolver::calcInvR,
                                  &ShallowWaterSolver::calcLambdaR,
                                  &ShallowWaterSolver::calcGradForR,
                                  tau, false);
        arma::vec inv(4);
        inv[0] = n.x;
        inv[1] = n.y;
        inv[2] = 2;
        inv[3] = R;
        invs.push_back(inv);
    }

    double lambdaQ = calcLambdaQ(h, u, n);
    if (lambdaQ > 0) {
        double Q = extrapolateInv(cell, node, n,
                                  &ShallowWaterSolver::calcInvQ,
                                  &ShallowWaterSolver::calcLambdaQ,
                                  &ShallowWaterSolver::calcGradForQ,
                                  tau, false);
        arma::vec inv(4);
        inv[0] = n.x;
        inv[1] = n.y;
        inv[2] = -2;
        inv[3] = Q;
        invs.push_back(inv);
    }

    double lambdaS = calcLambdaS(h, u, n);
    double S = extrapolateInv(cell, node, n,
                              &ShallowWaterSolver::calcInvS,
                              &ShallowWaterSolver::calcLambdaS,
                              &ShallowWaterSolver::calcGradForS,
                              tau, false);
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

std::vector<arma::vec> ShallowWaterSolver::getInvFromCellIntr(Node *node, Cell *cell, double avgH,
                                                              Vector avgU, double tau) {
    std::vector<arma::vec> invs;

    long centerNodeID = cell->centerNodeID;
    double h = mesh->s1[centerNodeID][0];
    Vector u = mesh->v1[centerNodeID][0];
    Vector n = cell->nodeToTransferVector[node->ID];

    double lambdaR = calcLambdaR(h, u, n);
    if (lambdaR > 0) {
        double R = interpolateInv(cell, node, avgH, avgU, n, &ShallowWaterSolver::calcInvR,
                                  &ShallowWaterSolver::calcLambdaR,
                                  &ShallowWaterSolver::calcGradForR, tau);
        arma::vec inv(4);
        inv[0] = n.x;
        inv[1] = n.y;
        inv[2] = 2;
        inv[3] = R;
        invs.push_back(inv);
    }

    double lambdaQ = calcLambdaQ(h, u, n);
    if (lambdaQ > 0) {
        double Q = interpolateInv(cell, node, avgH, avgU, n, &ShallowWaterSolver::calcInvQ,
                                  &ShallowWaterSolver::calcLambdaQ,
                                  &ShallowWaterSolver::calcGradForQ, tau);
        arma::vec inv(4);
        inv[0] = n.x;
        inv[1] = n.y;
        inv[2] = -2;
        inv[3] = Q;
        invs.push_back(inv);
    }

    double lambdaS = calcLambdaS(h, u, n);
    double S =
        interpolateInv(cell, node, avgH, avgU, n, &ShallowWaterSolver::calcInvS,
                       &ShallowWaterSolver::calcLambdaS, &ShallowWaterSolver::calcGradForS, tau);
    arma::vec inv(5);
    inv[0] = -n.y;
    inv[1] = n.x;
    inv[2] = 0;
    inv[3] = S;
    inv[4] = lambdaS;
    invs.push_back(inv);
    return invs;
}

double ShallowWaterSolver::extrapolateInv(
    Edge *edge, Node *node, Vector n,
    std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcInv,
    std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda,
    std::function<Vector(ShallowWaterSolver &, long, Vector)> calcGradForInv, double tau,
    bool needMonotize) {
    Data *data = &node->data;
    Data *oppositeData = &mesh->nodes[edge->getAnotherEndNode(node->ID)].data;

    Data *innerNode1 = &mesh->nodes[edge->getNearInnerNode(node->ID)].data;
    Data *innerNode2 = &mesh->nodes[edge->getFarInnerNode(node->ID)].data;
    double h1 = (innerNode1->getS0(0) + innerNode2->getS0(0) + innerNode1->getS2(0) +
                 innerNode2->getS2(0)) /
                4;
    Vector u1 = (innerNode1->getV0(0) + innerNode2->getV0(0) + innerNode1->getV2(0) +
                 innerNode2->getV2(0)) /
                4;

    // Cell *cell1 = &mesh->cells[edge->cellIDs[0]];
    // Cell *cell2 = &mesh->cells[edge->cellIDs[1]];
    // Data *centerNode1 = &mesh->nodes[cell1->ID].data;
    // Data *centerNode2 = &mesh->nodes[cell2->ID].data;
    // double h1 = (centerNode1->s1[0] + centerNode2->s1[0]) / 2;
    // Vector u1 = (centerNode1->v1[0] + centerNode2->v1[0]) / 2;

    double invCenter1 = calcInv(*this, h1, u1, n);
    double invOpposite0 = calcInv(*this, oppositeData->getS0(0), oppositeData->getV0(0), n);
    double coefPan = 0.0;
    double inv2 = (2 * invCenter1 - (1 - coefPan) * invOpposite0) / (1 + coefPan);
    // double inv2 = 2 * invCenter1 - invOpposite0;

    // double invInnerNode10 = calcInv(*this, innerNode1->s0[0], innerNode1->v0[0], n);
    // double invInnerNode20 = calcInv(*this, innerNode2->s0[0], innerNode2->v0[0], n);

    // double invInnerNode12 = calcInv(*this, innerNode1->s2[0], innerNode1->v2[0], n);
    // double invInnerNode22 = calcInv(*this, innerNode2->s2[0], innerNode2->v2[0], n);
    // double inv2 = 2 * invInnerNode12 - invInnerNode22;
    if (needMonotize) {
        double inv0 = calcInv(*this, data->getS0(0), data->getV0(0), n);
        double h0 = (innerNode1->getS0(0) + innerNode2->getS0(0)) / 2;
        Vector u0 = (innerNode1->getV0(0) + innerNode2->getV0(0)) / 2;
        // double h0 = (centerNode1->s0[0] + centerNode2->s0[0]) / 2;
        // Vector u0 = (centerNode1->v0[0] + centerNode2->v0[0]) / 2;
        double invCenter0 = calcInv(*this, h0, u0, n);

        double lambda = calcLambda(*this, h1, u1, n);
        Vector grad1 = calcGradForInv(*this, edge->cellIDs[0], n);
        Vector grad2 = calcGradForInv(*this, edge->cellIDs[1], n);
        Vector grad = (grad1 + grad2) / 2;

        inv2 =
            monotize(inv2, inv0, invCenter0, invOpposite0, invCenter1, lambda, tau, n * grad, true);

        // double lambda = calcLambda(*this, h1, u1, n);
        // double h = (innerNode1->coords - innerNode2->coords).length();
        // double dirDerivative = (invInnerNode10 - invInnerNode20 + invInnerNode12 -
        // invInnerNode22) / h;

        // inv2 = monotize(inv2, inv0, invCenter0, invOpposite0, invCenter1, lambda,
        //                 tau, dirDerivative, true);
    }

    return inv2;
}

std::vector<arma::vec> ShallowWaterSolver::getInvFromEdge(Node *node, Edge *edge, double avgH,
                                                          Vector avgU, double tau) {
    std::vector<arma::vec> invs;

    Data *innerNode1 = &mesh->nodes[edge->innerNodeIDs[0]].data;
    Data *innerNode2 = &mesh->nodes[edge->innerNodeIDs[1]].data;

    double h = (innerNode1->getS0(0) + innerNode2->getS0(0) + innerNode1->getS2(0) +
                innerNode2->getS2(0)) /
               4;
    Vector u = (innerNode1->getV0(0) + innerNode2->getV0(0) + innerNode1->getV2(0) +
                innerNode2->getV2(0)) /
               4;

    // Cell *cell1 = &mesh->cells[edge->cellIDs[0]];
    // Cell *cell2 = &mesh->cells[edge->cellIDs[1]];
    // Data *centerNode1 = &mesh->nodes[cell1->ID].data;
    // Data *centerNode2 = &mesh->nodes[cell2->ID].data;
    // double h = (centerNode1->s1[0] + centerNode2->s1[0]) / 2;
    // Vector u = (centerNode1->v1[0] + centerNode2->v1[0]) / 2;

    Node *oppositeNode = &mesh->nodes[edge->getAnotherEndNode(node->ID)];
    Vector n = node->data.coords - oppositeNode->data.coords;
    n = n / n.length();

    double lambdaR = calcLambdaR(h, u, n);
    if (lambdaR > 0) {
        double R = extrapolateInv(edge, node, n, &ShallowWaterSolver::calcInvR,
                                  &ShallowWaterSolver::calcLambdaR,
                                  &ShallowWaterSolver::calcGradForR, tau, true);
        arma::vec inv(4);
        inv[0] = n.x;
        inv[1] = n.y;
        inv[2] = 2;
        inv[3] = R;
        invs.push_back(inv);
    }

    double lambdaQ = calcLambdaQ(h, u, n);
    if (lambdaQ > 0) {
        double Q = extrapolateInv(edge, node, n, &ShallowWaterSolver::calcInvQ,
                                  &ShallowWaterSolver::calcLambdaQ,
                                  &ShallowWaterSolver::calcGradForQ, tau, true);
        arma::vec inv(4);
        inv[0] = n.x;
        inv[1] = n.y;
        inv[2] = -2;
        inv[3] = Q;
        invs.push_back(inv);
    }

    double lambdaS = calcLambdaS(h, u, n);
    double S = extrapolateInv(edge, node, n, &ShallowWaterSolver::calcInvS,
                              &ShallowWaterSolver::calcLambdaS, &ShallowWaterSolver::calcGradForS,
                              tau, true);
    arma::vec inv(5);
    inv[0] = -n.y;
    inv[1] = n.x;
    inv[2] = 0;
    inv[3] = S;
    inv[4] = lambdaS;
    invs.push_back(inv);
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
            Vector n = cell->nodeToTransferVector[nodeID];
            invs[1][3] = calcInvS(0, mesh->v1[centerNodeID][0], n);
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
        //     invs[1][3] = calcInvS(0, mesh->v1[centerNodeID][0], n);
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

Cell *ShallowWaterSolver::chooseCell(
    Node *node, double avgH, Vector avgU,
    std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda) {
    Cell *cellWithMaxLambda = nullptr;
    double maxLambda = -std::numeric_limits<double>::max();
    for (unsigned long cellID : node->cellIDs) {
        Cell *cell = &mesh->cells[cellID];
        Vector n = cell->nodeToTransferVector[node->ID];
        double lambda = calcLambda(*this, avgH, avgU, n);
        if (lambda >= maxLambda) {
            maxLambda = lambda;
            cellWithMaxLambda = cell;
        }
    }
    assert(cellWithMaxLambda != nullptr);

    return cellWithMaxLambda;
}

Cell *ShallowWaterSolver::chooseOppositeCell(Node *node, Cell *cell) {
    Cell *oppositeCell = nullptr;
    Vector n = cell->nodeToTransferVector[node->ID];
    Vector oppositeN = n;
    for (long cellID : node->cellIDs) {
        if (cellID == cell->ID) {
            continue;
        }

        Cell *tmpCell = &mesh->cells[cellID];
        Vector tmpN = tmpCell->nodeToTransferVector[node->ID];
        if (n * tmpN < n * oppositeN) {
            oppositeCell = tmpCell;
            oppositeN = tmpN;
        }
    }
    assert(oppositeCell != nullptr);

    return oppositeCell;
}

Edge *ShallowWaterSolver::chooseEdge(
    Node *node, double avgH, Vector avgU,
    std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda) {
    Edge *edgeWithMaxLambda = nullptr;
    double maxLambda = -std::numeric_limits<double>::max();
    for (unsigned long edgeID : node->edgeIDs) {
        Edge *edge = &mesh->edges[edgeID];

        Node *oppositeNode = &mesh->nodes[edge->getAnotherEndNode(node->ID)];
        Vector n = node->data.coords - oppositeNode->data.coords;
        n = n / n.length();

        double lambda = calcLambda(*this, avgH, avgU, n);
        if (lambda >= maxLambda) {
            maxLambda = lambda;
            edgeWithMaxLambda = edge;
        }
    }
    assert(edgeWithMaxLambda != nullptr);

    return edgeWithMaxLambda;
}

Edge *ShallowWaterSolver::chooseOppositeEdge(Node *node, Edge *edge) {
    Edge *oppositeEdge = nullptr;

    Node *oppositeNode = &mesh->nodes[edge->getAnotherEndNode(node->ID)];
    Vector n = node->data.coords - oppositeNode->data.coords;
    n = n / n.length();
    Vector oppositeN = n;
    for (long edgeID : node->edgeIDs) {
        if (edgeID == edge->ID) {
            continue;
        }

        Edge *tmpEdge = &mesh->edges[edgeID];

        Node *tmpOppositeNode = &mesh->nodes[tmpEdge->getAnotherEndNode(node->ID)];
        Vector tmpN = node->data.coords - tmpOppositeNode->data.coords;
        tmpN = tmpN / tmpN.length();
        if (n * tmpN < n * oppositeN) {
            oppositeEdge = tmpEdge;
            oppositeN = tmpN;
        }
    }
    assert(oppositeEdge != nullptr);

    return oppositeEdge;
}

void ShallowWaterSolver::processPhase2BoundNode(Node *node, double tau) {
    double avgH = 0;
    Vector avgU = Vector(0, 0);
    for (unsigned long edgeID : node->edgeIDs) {
        Edge *edge = &mesh->edges[edgeID];
        long neighborInnerNodeID = edge->getNearInnerNode(node->ID);
        avgH = avgH + mesh->s2[neighborInnerNodeID][0];
        avgU = avgU + mesh->v2[neighborInnerNodeID][0];
    }
    avgH = avgH / node->edgeIDs.size();
    avgU = avgU / node->edgeIDs.size();

    Cell *cell = chooseCell(node, avgH, avgU, &ShallowWaterSolver::calcLambdaR);

    std::vector<arma::vec> invs = getInvFromCellIntr(node, cell, avgH, avgU, tau);
    assert(invs.size() == 2);

    // --- no-slip condition ---
    // mesh->v2[node->ID][0].x = 0;
    // mesh->v2[node->ID][0].y = 0;
    // mesh->s2[node->ID][0] = pow(invs[0][3] / 2, 2) / g;
    // --- ---

    // --- slip condition ---
    std::vector<Vector> edgeNs;
    for (unsigned long edgeID : node->edgeIDs) {
        Edge *edge = &mesh->edges[edgeID];
        if (edge->boundEdge) {
            edgeNs.push_back(edge->normal);
        }
    }
    assert(edgeNs.size() == 2);
    double cos = edgeNs[0] * edgeNs[1];
    if (cos >= 1 - 0.1e-9) {
        Vector edgeN = edgeNs[0];
        arma::vec inv(4);
        inv[0] = edgeN.x;
        inv[1] = edgeN.y;
        inv[2] = 0;
        inv[3] = 0;
        invs.push_back(inv);

        if (invs[1][4] <= 0) {
            Vector n = cell->nodeToTransferVector[node->ID];
            invs[1][3] = calcInvS(0, mesh->v1[cell->centerNodeID][0], n);
        }

        arma::vec initialValues = convertInvToInitialVariables(invs);

        mesh->v2[node->ID][0].x = initialValues[0];
        mesh->v2[node->ID][0].y = initialValues[1];
        mesh->s2[node->ID][0] = initialValues[2];
    } else {
        mesh->v2[node->ID][0].x = 0;
        mesh->v2[node->ID][0].y = 0;
        // mesh->s2[node->ID][0] = mesh->s1[cell->centerNodeID][0];
        double avgH = 0;
        for (unsigned long edgeID : node->edgeIDs) {
            Edge *edge = &mesh->edges[edgeID];
            avgH += mesh->s0[edge->getNearInnerNode(node->ID)][0];
        }
        mesh->s2[node->ID][0] = avgH / node->edgeIDs.size();
    }
    // --- ---

    // --- flow condition ---
    // std::vector<Vector> edgeNs;
    // for (unsigned long edgeID : node->edgeIDs) {
    //     Edge *edge = &mesh->edges[edgeID];
    //     if (edge->boundEdge) {
    //         edgeNs.push_back(edge->normal);
    //     }
    // }
    // if (edgeNs.size() != 2) {
    //     std::cout << edgeNs.size() << std::endl;
    // }
    // assert(edgeNs.size() == 2);
    // double cos = edgeNs[0] * edgeNs[1];
    // if (cos >= 1 - 0.1e-9) {
    //     Vector edgeN = edgeNs[0];
    //     arma::vec inv(4);
    //     inv[0] = -edgeN.x;
    //     inv[1] = -edgeN.y;
    //     inv[2] = 2;
    //     inv[3] = 2 * sqrt(g);
    //     invs.push_back(inv);

    //     if (invs[1][4] <= 0) {
    //         Vector n = cell->nodeToTransferVector[node->ID];
    //         invs[1][3] = calcInvS(0, mesh->v1[cell->centerNodeID][0], n);
    //     }

    //     arma::vec initialValues = convertInvToInitialVariables(invs);

    //     mesh->v2[node->ID][0].x = initialValues[0];
    //     mesh->v2[node->ID][0].y = initialValues[1];
    //     mesh->s2[node->ID][0] = initialValues[2];
    // } else {
    //     double avgUx = 0;
    //     double avgUy = 0;
    //     double avgH = 0;
    //     for (unsigned long edgeID : node->edgeIDs) {
    //         Edge *edge = &mesh->edges[edgeID];
    //         long nearestInnerNodeId = edge->getNearInnerNode(node->ID);
    //         avgH += mesh->s0[nearestInnerNodeId][0];
    //         avgUx += mesh->v0[nearestInnerNodeId][0].x;
    //         avgUy += mesh->v0[nearestInnerNodeId][0].y;
    //     }
    //     mesh->v2[node->ID][0].x = avgUx / node->edgeIDs.size();
    //     mesh->v2[node->ID][0].y = avgUy / node->edgeIDs.size();
    //     mesh->s2[node->ID][0] = avgH / node->edgeIDs.size();
    // }
    // --- ---
}

void ShallowWaterSolver::processPhase2InnerNode(Node *node, double tau) {
    // --- calculate new values base on avg values in nearest neig ---
    // double avgH = 0;
    // Vector avgU = Vector(0, 0);
    // for (unsigned long edgeID : node->edgeIDs) {
    //     Edge *edge = &mesh->edges[edgeID];
    //     long neighborInnerNodeID = edge->getNearInnerNode(node->ID);
    //     Data *data = &mesh->nodes[neighborInnerNodeID].data;
    //     avgH = avgH + data->s2[0];
    //     avgU = avgU + data->v2[0];
    // }
    // avgH = avgH / node->edgeIDs.size();
    // avgU = avgU / node->edgeIDs.size();

    // Data *data = &node->data;
    // data->s2[0] = avgH;
    // data->v2[0] = avgU;
    // -------

    // --- calculate new values by interpolation in cell ---
    double avgH = 0;
    Vector avgU = Vector(0, 0);
    for (unsigned long edgeID : node->edgeIDs) {
        Edge *edge = &mesh->edges[edgeID];
        long neighborInnerNodeID = edge->getNearInnerNode(node->ID);
        avgH = avgH + mesh->s2[neighborInnerNodeID][0];
        avgU = avgU + mesh->v2[neighborInnerNodeID][0];
    }
    avgH = avgH / node->edgeIDs.size();
    avgU = avgU / node->edgeIDs.size();

    Cell *cell1 = chooseCell(node, avgH, avgU, &ShallowWaterSolver::calcLambdaR);
    Cell *cell2 = chooseOppositeCell(node, cell1);

    std::vector<arma::vec> invs1 = getInvFromCellIntr(node, cell1, avgH, avgU, tau);
    std::vector<arma::vec> invs2 = getInvFromCellIntr(node, cell2, avgH, avgU, tau);

    assert(invs1.size() == 2);
    assert(invs2.size() == 2);

    std::vector<arma::vec> invs;
    invs.push_back(invs1[0]);
    invs.push_back(invs2[0]);
    if (invs1[1][4] > invs2[1][4]) {
        invs.push_back(invs1[1]);
    } else {
        invs.push_back(invs2[1]);
    }
    arma::vec initialValues = convertInvToInitialVariables(invs);

    mesh->v2[node->ID][0].x = initialValues[0];
    mesh->v2[node->ID][0].y = initialValues[1];
    mesh->s2[node->ID][0] = initialValues[2];
    // -------

    // --- calculate new values by using values from edges ---
    // double avgH = 0;
    // Vector avgU = Vector(0, 0);
    // for (unsigned long edgeID : node->edgeIDs) {
    //     Edge *edge = &mesh->edges[edgeID];
    //     long neighborInnerNodeID = edge->getNearInnerNode(node->ID);
    //     Data *data = &mesh->nodes[neighborInnerNodeID].data;
    //     avgH = avgH + data->s2[0];
    //     avgU = avgU + data->v2[0];
    // }
    // avgH = avgH / node->edgeIDs.size();
    // avgU = avgU / node->edgeIDs.size();

    // Edge *edge1 = chooseEdge(node, avgH, avgU, &ShallowWaterSolver::calcLambdaR);
    // Edge *edge2 = chooseOppositeEdge(node, edge1);

    // std::vector<arma::vec> invs1 = getInvFromEdge(node, edge1, avgH, avgU, tau);
    // std::vector<arma::vec> invs2 = getInvFromEdge(node, edge2, avgH, avgU, tau);

    // assert(invs1.size() == 2);
    // assert(invs2.size() == 2);

    // std::vector<arma::vec> invs;
    // invs.push_back(invs1[0]);
    // invs.push_back(invs2[0]);
    // if (invs1[1][4] > invs2[1][4]) {
    //     invs.push_back(invs1[1]);
    // } else {
    //     invs.push_back(invs2[1]);
    // }
    // arma::vec initialValues = convertInvToInitialVariables(invs);

    // Data *data = &node->data;
    // data->v2[0].x = initialValues[0];
    // data->v2[0].y = initialValues[1];
    // data->s2[0] = initialValues[2];
    // -------
}

void ShallowWaterSolver::processPhase2(double tau) {
    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];

        Vector sqrtGH = calcGrad(cell, &ShallowWaterSolver::get2SqrtGH);
        Vector ux = calcGrad(cell, &ShallowWaterSolver::getUx);
        Vector uy = calcGrad(cell, &ShallowWaterSolver::getUy);

        mesh->cellIDToGrads[cell->ID][0].x = sqrtGH.x;
        mesh->cellIDToGrads[cell->ID][0].y = sqrtGH.y;
        mesh->cellIDToGrads[cell->ID][1].x = ux.x;
        mesh->cellIDToGrads[cell->ID][1].y = ux.y;
        mesh->cellIDToGrads[cell->ID][2].x = uy.x;
        mesh->cellIDToGrads[cell->ID][2].y = uy.y;
    }

    for (unsigned long i = 0; i < mesh->edges.size(); i++) {
        Edge *edge = &mesh->edges[i];
        if (edge->boundEdge) {
            processPhase2BoundEdge(edge, tau);
        } else {
            processPhase2InnerEdge(edge, tau);
        }
    }

    for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
        Node *node = &mesh->nodes[i];
        if (!node->isApex || !node->used) {
            continue;
        }
        if (node->boundNode) {
            processPhase2BoundNode(node, tau);
        } else {
            processPhase2InnerNode(node, tau);
        }
    }
}

void ShallowWaterSolver::processPhase3(double tau) {
    for (unsigned long i = 0; i < mesh->edges.size(); i++) {
        Edge *edge = &mesh->edges[i];

        std::vector<double> div1Edge;
        std::vector<double> div2Edge;
        std::vector<double> div3Edge;

        Vector n = edge->normal;

        for (unsigned long usedNodeID : edge->getUsedNodes(*mesh)) {
            double h = mesh->s2[usedNodeID][0];
            Vector u = mesh->v2[usedNodeID][0];

            div1Edge.push_back(calcDiv1(h, u, n));
            div2Edge.push_back(calcDiv2(h, u, n));
            div3Edge.push_back(calcDiv3(h, u, n));
        }

        mesh->edgeIDToDivs[edge->ID][0] = calcIntegral(div1Edge, edge->length);
        mesh->edgeIDToDivs[edge->ID][1] = calcIntegral(div2Edge, edge->length);
        mesh->edgeIDToDivs[edge->ID][2] = calcIntegral(div3Edge, edge->length);
    }

    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];

        double div1 = 0;
        double div2 = 0;
        double div3 = 0;
        for (long edgeID : cell->edgeIDs) {
            int coef = cell->edgeToNormalDir[edgeID];

            div1 += coef * mesh->edgeIDToDivs[edgeID][0];
            div2 += coef * mesh->edgeIDToDivs[edgeID][1];
            div3 += coef * mesh->edgeIDToDivs[edgeID][2];
        }
        div1 /= cell->volume;
        div2 /= cell->volume;
        div3 /= cell->volume;

        long centerNodeID = cell->centerNodeID;

        double h = mesh->s1[centerNodeID][0];
        Vector u = mesh->v1[centerNodeID][0];
        double newH = h - tau * div1 / 2;
        double newUX = (h * u.x - tau * div2 / 2) / newH;
        double newUY = (h * u.y - tau * div3 / 2) / newH;

        mesh->s2[centerNodeID][0] = newH;
        mesh->v2[centerNodeID][0].x = newUX;
        mesh->v2[centerNodeID][0].y = newUY;
    }
}

void ShallowWaterSolver::prepareNextStep() {
    for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
        for (unsigned long j = 0; j < mesh->s0[i].size(); j++) {
            mesh->s0[i][j] = mesh->s2[i][j];
        }
        for (unsigned long j = 0; j < mesh->v0[i].size(); j++) {
            mesh->v0[i][j] = mesh->v2[i][j];
        }
    }
}