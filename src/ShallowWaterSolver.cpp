#include "ShallowWaterSolver.hpp"

#include <armadillo>
#include <cassert>
#include <iomanip>
#include <limits>
#include <cmath>

double ShallowWaterSolver::calcTau() {
    double tau = std::numeric_limits<double>::max();
    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Data *cellData = &mesh->nodes[cell->centerNodeID].data;
        for (unsigned long edgeID : cell->edgeIDs) {
            for (unsigned long nodeID : mesh->edges[edgeID].innerNodeIDs) {
                Data *data = &mesh->nodes[nodeID].data;
                
                double h = 2 * (data->coords - cellData->coords).length();
                double lambdaR = calcLambdaR(data->s0[0], data->v0[0], 
                                            cell->nodeToTransferVector[nodeID]);
                double nodeTauR = h / std::abs(lambdaR);
                if (tau > nodeTauR) {
                    tau = nodeTauR;
                }
                double lambdaQ = calcLambdaQ(data->s0[0], data->v0[0], 
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
    return u * n * h;
}

double ShallowWaterSolver::calcDiv2(double h, Vector u, Vector n) {
    Vector F = Vector(pow(u.x, 2) + g * h / 2, u.x * u.y);
    return F * n * h;
}

double ShallowWaterSolver::calcDiv3(double h, Vector u, Vector n) {
    Vector F = Vector(u.x * u.y, pow(u.y, 2) + g * h / 2);
    return F * n * h;
}

double calcIntegral(std::vector<double> values, double length) {
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
            Data *data = &mesh->nodes[usedNodeID].data;
            double h = data->s0[0];
            Vector u = data->v0[0];

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

        Data *cellData = &mesh->nodes[cell->centerNodeID].data;

        double h = cellData->s0[0];
        Vector u = cellData->v0[0];
        double newH = h - tau * div1 / 2;
        double newUX = (h * u.x - tau * div2 / 2) / newH;
        double newUY = (h * u.y - tau * div3 / 2) / newH;

        cellData->s1[0] = newH;
        cellData->v1[0].x = newUX;
        cellData->v1[0].y = newUY;
    }
}

double ShallowWaterSolver::calcInvR(double h, Vector u, double G, Vector n) {
    return u * n + 2 * sqrt(g * h); 
}

double ShallowWaterSolver::calcInvQ(double h, Vector u, double G, Vector n) {
    return u * n - 2 * sqrt(g * h); 
}

double ShallowWaterSolver::calcInvS(double h, Vector u, double G, Vector n) { 
    return -u.x * n.y + u.y * n.x; 
}

double ShallowWaterSolver::calcLambdaR(double h, Vector u, Vector n) { 
    return u * n + sqrt(g * h); 
}

double ShallowWaterSolver::calcLambdaQ(double h, Vector u, Vector n) { 
    return u * n - sqrt(g * h); 
}

double ShallowWaterSolver::calcLambdaS(double h, Vector u, Vector n) { 
    return u * n; 
}

double monotize(double inv2, double inv0, double invCenter0, double invOpposite0,
                double invCenter1, double lambda, double tau, double dirDerivative, 
                bool withQ) {

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

double ShallowWaterSolver::get2SqrtGH(Data *data) {
    return 2 * sqrt(g * data->s0[0]);
}

double ShallowWaterSolver::getUx(Data *data) {
    return data->v0[0].x;
}

double ShallowWaterSolver::getUy(Data *data) {
    return data->v0[0].y;
}

Vector ShallowWaterSolver::calcGrad(
    Cell *cell, std::function<double(ShallowWaterSolver &, Data*)> getVar) {

    Vector grad = Vector(0, 0);

    for (long edgeID : cell->edgeIDs) {
        Edge *edge = &mesh->edges[edgeID];
        std::vector<double> variables;
        for (unsigned long usedNodeID : edge->getUsedNodes(*mesh)) {
            Data *data = &mesh->nodes[usedNodeID].data;
            variables.push_back(getVar(*this, data));
        }
        Vector edgeN = edge->normal * cell->edgeToNormalDir[edgeID] ;
        grad = grad +  edgeN * calcIntegral(variables, edge->length);
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
    Cell *cell, Node *node, Vector n, double G,
    std::function<double(ShallowWaterSolver &, double, Vector, double, Vector)> calcInv,
    std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda, 
    std::function<Vector(ShallowWaterSolver &, long, Vector)> calcGradForInv,
    double tau, bool needMonotize) {

    assert(cell->nodeIDToOppositeNodeID[node->ID] != -1);

    Data *data = &node->data;
    Data *centerData = &mesh->nodes[cell->centerNodeID].data;
    Data *oppositeData = &mesh->nodes[cell->nodeIDToOppositeNodeID[node->ID]].data;

    double invCenter1 = calcInv(*this, centerData->s1[0], centerData->v1[0], G, n);
    double invOpposite0 = calcInv(*this, oppositeData->s0[0], oppositeData->v0[0], G, n);
    double inv2 = 2 * invCenter1 - invOpposite0;
    if (needMonotize) {
        double inv0 = calcInv(*this, data->s0[0], data->v0[0], G, n);
        double invCenter0 = calcInv(*this, centerData->s0[0], centerData->v0[0], G, n);

        double lambda = calcLambda(*this, centerData->s1[0], centerData->v1[0], n);
        Vector grad = calcGradForInv(*this, cell->ID, n);
        inv2 = monotize(inv2, inv0, invCenter0, invOpposite0, invCenter1, lambda, 
                        tau, n * grad, true);
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
    Cell *cell, Node *node, double avgH, Vector avgU, Vector n, double G,
    std::function<double(ShallowWaterSolver &, double, Vector, double, Vector)> calcInv,
    std::function<double(ShallowWaterSolver &, double, Vector, Vector)> calcLambda,
    std::function<Vector(ShallowWaterSolver &, long, Vector)> calcGradForInv, 
    double tau) {

    arma::vec invs(10);
    int pos = 0;
    for (long nodeID : cell->nodeIDs) {
        Data *data = &mesh->nodes[nodeID].data;
        invs[pos] = calcInv(*this, data->s0[0], data->v0[0], G, n);
        pos++;
    }
    for (long edgeID : cell->edgeIDs) {
        for (long nodeID : mesh->edges[edgeID].innerNodeIDs) {
            Data *data = &mesh->nodes[nodeID].data;
            invs[pos] = calcInv(*this, data->s0[0], data->v0[0], G, n);
            pos++;
        }
    }
    Data *centerData = &mesh->nodes[cell->centerNodeID].data;
    invs[pos] = calcInv(*this, centerData->s0[0], centerData->v0[0], G, n);
    arma::vec inv = cell->interpolationMat * invs;

    Data *data = &mesh->nodes[node->ID].data;
    double coef = -calcLambda(*this, avgH, avgU, n) * tau;
    Vector intersectCoords = n * coef + data->coords;
    double intersectX = intersectCoords.x;
    double intersectY = intersectCoords.y;

    double inv2 = interpolate(intersectX, intersectY, inv);

    double invCenter0 = calcInv(*this, centerData->s0[0], centerData->v0[0], G, n);
    double invCenter1 = calcInv(*this, centerData->s1[0], centerData->v1[0], G, n);
    double lambda = calcLambda(*this, centerData->s1[0], centerData->v1[0], n);

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

    assert(arma::rank(A) == 3);
    arma::vec initialVariables = arma::solve(A, b);
    initialVariables[2] = pow(initialVariables[2], 2) / g;
    return initialVariables;
}

std::vector<arma::vec> ShallowWaterSolver::getInvFromCellExtr(Node *node, Cell *cell, double tau) {
    std::vector<arma::vec> invs;

    Data *centerData = &mesh->nodes[cell->centerNodeID].data;
    double h = centerData->s1[0];
    Vector u = centerData->v1[0];
    Vector n = cell->nodeToTransferVector[node->ID];
    
    double G = sqrt(g / h);

    double lambdaR = calcLambdaR(h, u, n);
    if (lambdaR > 0) {
        double R = extrapolateInv(cell, node, n, G, 
                                  &ShallowWaterSolver::calcInvR,
                                  &ShallowWaterSolver::calcLambdaR, 
                                  &ShallowWaterSolver::calcGradForR,
                                  tau, true);
        arma::vec inv(4);
        inv[0] = n.x;
        inv[1] = n.y;
        inv[2] = 2;
        inv[3] = R;
        invs.push_back(inv);
    }

    double lambdaQ = calcLambdaQ(h, u, n);
    if (lambdaQ > 0) {
        double Q = extrapolateInv(cell, node, n, G, 
                                  &ShallowWaterSolver::calcInvQ,
                                  &ShallowWaterSolver::calcLambdaQ, 
                                  &ShallowWaterSolver::calcGradForQ,
                                  tau, true);
        arma::vec inv(4);
        inv[0] = n.x;
        inv[1] = n.y;
        inv[2] = -2;
        inv[3] = Q;
        invs.push_back(inv);
    }

    double lambdaS = calcLambdaS(h, u, n);
    double S = extrapolateInv(cell, node, n, G,
                              &ShallowWaterSolver::calcInvS,
                              &ShallowWaterSolver::calcLambdaS, 
                              &ShallowWaterSolver::calcGradForS,
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

std::vector<arma::vec> ShallowWaterSolver::getInvFromCellIntr(Node *node, Cell *cell,
                                                              double avgH, Vector avgU, 
                                                              double tau) {
    std::vector<arma::vec> invs;

    Data *centerData = &mesh->nodes[cell->centerNodeID].data;
    double h = centerData->s1[0];
    Vector u = centerData->v1[0];
    Vector n = cell->nodeToTransferVector[node->ID];
    
    double G = sqrt(g / h);

    double lambdaR = calcLambdaR(h, u, n);
    if (lambdaR > 0) {
        double R = interpolateInv(cell, node, avgH, avgU, n, G, 
                                  &ShallowWaterSolver::calcInvR,
                                  &ShallowWaterSolver::calcLambdaR, 
                                  &ShallowWaterSolver::calcGradForR,
                                  tau);
        arma::vec inv(4);
        inv[0] = n.x;
        inv[1] = n.y;
        inv[2] = 2;
        inv[3] = R;
        invs.push_back(inv);
    }

    double lambdaQ = calcLambdaQ(h, u, n);
    if (lambdaQ > 0) {
        double Q = interpolateInv(cell, node, avgH, avgU, n, G, 
                                  &ShallowWaterSolver::calcInvQ,
                                  &ShallowWaterSolver::calcLambdaQ, 
                                  &ShallowWaterSolver::calcGradForQ,
                                  tau);
        arma::vec inv(4);
        inv[0] = n.x;
        inv[1] = n.y;
        inv[2] = -2;
        inv[3] = Q;
        invs.push_back(inv);
    }

    double lambdaS = calcLambdaS(h, u, n);
    double S = interpolateInv(cell, node, avgH, avgU, n, G, 
                              &ShallowWaterSolver::calcInvS,
                              &ShallowWaterSolver::calcLambdaS, 
                              &ShallowWaterSolver::calcGradForS,
                              tau);
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

    for (unsigned long nodeID : edge->innerNodeIDs) {
        Node *node = &mesh->nodes[nodeID];

        std::vector<arma::vec> invs = getInvFromCellExtr(node, cell, tau);
        assert(invs.size() == 2);

        Data *data = &node->data;
        data->v2[0].x = 0;
        data->v2[0].y = 0;
        data->s2[0] = pow(invs[0][3] / 2, 2) / g;

        // Vector n = cell->nodeToTransferVector[node->ID];
        // arma::vec inv(4);
        // inv[0] = -n.x;
        // inv[1] = -n.y;
        // inv[2] = 2;
        // inv[3] = 2 * sqrt(g * 1);
        // invs.push_back(inv);

        // if (invs[1][4] <= 0) {
        //     Vector n = cell->nodeToTransferVector[nodeID];
        //     invs[1][0] = n.y;
        //     invs[1][1] = -n.x;
        //     invs[1][2] = 0;
        //     invs[1][3] = 0;
        // }

        // if (invs[1][4] <= 0) {
        //     Data *cellData = &mesh->nodes[cell->centerNodeID].data;
        //     Vector n = cell->nodeToTransferVector[node->ID];
        //     invs[1][3] = calcInvS(cellData->s1[0], cellData->v1[0], 0, n);
        // }

        // arma::vec initialValues = convertInvToInitialVariables(invs);

        // Data *data = &node->data;
        // data->v2[0].x = initialValues[0];
        // data->v2[0].y = initialValues[1];
        // data->s2[0] = initialValues[2];
    }
}

void ShallowWaterSolver::processPhase2InnerEdge(Edge *edge, double tau) {
    Cell *cell1 = &mesh->cells[edge->cellIDs[0]];
    Cell *cell2 = &mesh->cells[edge->cellIDs[1]];

    for (unsigned long nodeID : edge->innerNodeIDs) {
        Node *node = &mesh->nodes[nodeID];

        std::vector<arma::vec> invs1 = getInvFromCellExtr(node, cell1, tau);
        std::vector<arma::vec> invs2 = getInvFromCellExtr(node, cell2, tau);
        std::vector<arma::vec> invs;

        assert(invs1.size() == 2);
        assert(invs2.size() == 2);

        invs.push_back(invs1[0]);
        invs.push_back(invs2[0]);
        if (invs1[1][4] > invs2[1][4]) {
            invs.push_back(invs1[1]);
        } else {
            invs.push_back(invs2[1]);
        }
        arma::vec initialValues = convertInvToInitialVariables(invs);

        Data *data = &node->data;
        data->v2[0].x = initialValues[0];
        data->v2[0].y = initialValues[1];
        data->s2[0] = initialValues[2];
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
    for (unsigned long cellID : node->cellIDs) {
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

void ShallowWaterSolver::processPhase2BoundNode(Node *node, double tau) {
    double avgH = 0;
    Vector avgU = Vector(0, 0);
    for (unsigned long edgeID : node->edgeIDs) {
        Edge *edge = &mesh->edges[edgeID];
        long neighborInnerNodeID = edge->getNearInnerNode(node->ID);
        Data *data = &mesh->nodes[neighborInnerNodeID].data;
        avgH = avgH + data->s2[0];
        avgU = avgU + data->v2[0];
    }
    avgH = avgH / node->edgeIDs.size();
    avgU = avgU / node->edgeIDs.size();

    Cell *cell = chooseCell(node, avgH, avgU, &ShallowWaterSolver::calcLambdaR);

    std::vector<arma::vec> invs = getInvFromCellIntr(node, cell, avgH, avgU, tau);
    assert(invs.size() == 2);

    Data *data = &node->data;
    data->v2[0].x = 0;
    data->v2[0].y = 0;
    data->s2[0] = pow(invs[0][3] / 2, 2) / g;
        
    // Vector n = cell->nodeToTransferVector[node->ID];
    // arma::vec inv(4);
    // inv[0] = -n.x;
    // inv[1] = -n.y;
    // inv[2] = 2;
    // inv[3] = 2 * sqrt(g * 1);
    // invs.push_back(inv);

    // if (invs[1][4] <= 0) {
    //     Vector n = cell->nodeToTransferVector[node->ID];
    //     invs[1][0] = n.y;
    //     invs[1][1] = -n.x;
    //     invs[1][2] = 0;
    //     invs[1][3] = 0;
    // }

    // if (invs[1][4] <= 0) {
    //     Data *cellData = &mesh->nodes[cell->centerNodeID].data;
    //     Vector n = cell->nodeToTransferVector[node->ID];
    //     invs[1][3] = calcInvS(cellData->s1[0], cellData->v1[0], 0, n);
    // }

    // arma::vec initialValues = convertInvToInitialVariables(invs);

    // Data *data = &node->data;
    // data->v2[0].x = initialValues[0];
    // data->v2[0].y = initialValues[1];
    // data->s2[0] = initialValues[2];
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
        Data *data = &mesh->nodes[neighborInnerNodeID].data;
        avgH = avgH + data->s2[0];
        avgU = avgU + data->v2[0];
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

    Data *data = &node->data;
    data->v2[0].x = initialValues[0];
    data->v2[0].y = initialValues[1];
    data->s2[0] = initialValues[2];
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
            Data *data = &mesh->nodes[usedNodeID].data;
            double h = data->s2[0];
            Vector u = data->v2[0];

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

        Data *cellData = &mesh->nodes[cell->centerNodeID].data;

        double h = cellData->s1[0];
        Vector u = cellData->v1[0];
        double newH = h - tau * div1 / 2;
        double newUX = (h * u.x - tau * div2 / 2) / newH;
        double newUY = (h * u.y - tau * div3 / 2) / newH;

        cellData->s2[0] = newH;
        cellData->v2[0].x = newUX;
        cellData->v2[0].y = newUY;
    }
}

void ShallowWaterSolver::prepareNextStep() {
    for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
        Data *data = &mesh->nodes[i].data;

        for (unsigned long j = 0; j < data->s0.size(); j++) {
            data->s0[j] = data->s2[j];
        }
        for (unsigned long j = 0; j < data->v0.size(); j++) {
            data->v0[j] = data->v2[j];
        }
    }
}