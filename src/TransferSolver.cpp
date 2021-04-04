// #include "TransferSolver.hpp"

// #include <armadillo>
// #include <cassert>
// #include <iomanip>
// #include <limits>

// double TransferSolver::calcTau() {
    // double tau = std::numeric_limits<double>::max();
    // for (unsigned long i = 0; i < mesh->cells.size(); i++) {
    //     Cell *cell = &mesh->cells[i];
    //     Node *cellNode = &mesh->nodes[cell->centerNodeID];
    //     for (unsigned long edgeID : cell->edgeIDs) {
    //         for (unsigned long nodeID : mesh->edges[edgeID].innerNodeIDs) {
    //             Node *node = &mesh->nodes[nodeID];
    //             double lambda = node->data.vector[0] * cell->nodeToTransferVector[nodeID];
    //             double h = 2 * (node->data.coords - cellNode->data.coords).length();
    //             double nodeTau = h / std::abs(lambda);
    //             if (tau > nodeTau) {
    //                 tau = nodeTau;
    //             }
    //         }
    //     }
    // }
    // assert(tau != std::numeric_limits<double>::max());
    // return cfl * tau;
// }

// double TransferSolver::calcDivOnEdge(Edge *edge, int phase) {
//     std::vector<long> usedNodeIDs = edge->usedNodeIDs;
//     double div = 0;
//     switch (usedNodeIDs.size()) {
//         case 1: {
//             Node *edgeNode = &mesh->nodes[usedNodeIDs[0]];

//             double phi = phase == 1 ? edgeNode->data.s0[0] : edgeNode->data.s2[0];
//             div = phi * (edgeNode->data.vector[0] * edge->normal) * edge->length;
//             break;
//         }
//         case 2: {
//             Node *edgeNode1 = &mesh->nodes[usedNodeIDs[0]];
//             Node *edgeNode2 = &mesh->nodes[usedNodeIDs[1]];

//             double phi1 = phase == 1 ? edgeNode1->data.s0[0] : edgeNode1->data.s2[0];
//             double phi2 = phase == 1 ? edgeNode2->data.s0[0] : edgeNode2->data.s2[0];
//             double phi = (phi1 + phi2) / 2;
//             div = phi * (edgeNode1->data.vector[0] * edge->normal) * edge->length;
//             break;
//         }
//         case 4: {
//             Node *edgeNode1 = &mesh->nodes[usedNodeIDs[0]];
//             Node *edgeNode2 = &mesh->nodes[usedNodeIDs[1]];
//             Node *edgeNode3 = &mesh->nodes[usedNodeIDs[2]];
//             Node *edgeNode4 = &mesh->nodes[usedNodeIDs[3]];

//             double phi1 = phase == 1 ? edgeNode1->data.s0[0] : edgeNode1->data.s2[0];
//             double phi2 = phase == 1 ? edgeNode2->data.s0[0] : edgeNode2->data.s2[0];
//             double phi3 = phase == 1 ? edgeNode3->data.s0[0] : edgeNode3->data.s2[0];
//             double phi4 = phase == 1 ? edgeNode4->data.s0[0] : edgeNode4->data.s2[0];

//             double div1 = phi1 * (edgeNode1->data.vector[0] * edge->normal);
//             double div2 = phi2 * (edgeNode2->data.vector[0] * edge->normal);
//             double div3 = phi3 * (edgeNode3->data.vector[0] * edge->normal);
//             double div4 = phi4 * (edgeNode4->data.vector[0] * edge->normal);

//             div = (div1 + 3 * div2 + 3 * div3 + div4) / 8 * edge->length;
//             break;
//         }
//         default:
//             assert(false);
//             break;
//     }

//     return div;
// }

// void TransferSolver::processPhase1(double tau) {
//     for (unsigned long i = 0; i < mesh->cells.size(); i++) {
//         Cell *cell = &mesh->cells[i];
//         Node *cellCenter = &mesh->nodes[cell->centerNodeID];
//         double div = 0;

//         for (long edgeID : cell->edgeIDs) {
//             div += calcDivOnEdge(&mesh->edges[edgeID], 1) * cell->edgeToNormalDir[edgeID];
//         }
//         div /= cell->volume;
//         cellCenter->data.s1[0] = cellCenter->data.s0[0] - tau * div / 2;
//     }
// }

// double TransferSolver::getNewInvariantValue(Data *data, Data *centerData, Data *oppositeData,
//                                             Vector transfer, double tau, bool monotize) {
//     // double phiOpposite0 = oppositeData->s0[0];
//     double phiCenter0 = centerData->s0[0];
//     double phi0 = data->s0[0];
//     double phiOpposite0 = 2 * phiCenter0 - phi0;

//     double phiCenter1 = centerData->s1[0];

//     double phi2 = 2 * phiCenter1 - phiOpposite0;

//     if (monotize) {
//         double h = (oppositeData->coords - data->coords).length();
//         double Q = (phiCenter1 - phiCenter0) / (tau / 2) +
//                    (centerData->vector[0] * transfer) * (phi0 - phiOpposite0) / h;
//         double min = std::min(std::min(phi0, phiOpposite0), phiCenter0) + tau * Q;
//         double max = std::max(std::max(phi0, phiOpposite0), phiCenter0) + tau * Q;

//         phi2 = std::max(phi2, min);
//         phi2 = std::min(phi2, max);
//     }

//     return phi2;
// }

// void TransferSolver::processPhase2BoundEdge(Edge *edge, double tau) {
//     Cell *cell = &mesh->cells[edge->cellIDs[0]];
//     Data *centerData = &mesh->nodes[cell->centerNodeID].data;

//     for (unsigned long nodeID : edge->innerNodeIDs) {
//         Node *node = &mesh->nodes[nodeID];
//         Data *data = &node->data;
//         if (centerData->vector[0] * cell->nodeToTransferVector[nodeID] < 0) {
//             data->s2[0] = 0;
//             continue;
//         }

//         assert(cell->nodeIDToOppositeNodeID[nodeID] != -1);
//         Data *oppositeData = &mesh->nodes[cell->nodeIDToOppositeNodeID[nodeID]].data;

//         data->s2[0] = getNewInvariantValue(data, centerData, oppositeData,
//                                                  cell->nodeToTransferVector[nodeID], tau, true);
//     }
// }

// void TransferSolver::processPhase2InnerEdge(Edge *edge, double tau) {
//     Cell *cell1 = &mesh->cells[edge->cellIDs[0]];
//     Cell *cell2 = &mesh->cells[edge->cellIDs[1]];

//     Data *centerData1 = &mesh->nodes[cell1->centerNodeID].data;
//     Data *centerData2 = &mesh->nodes[cell2->centerNodeID].data;

//     Vector normal(edge->normal * cell1->edgeToNormalDir[edge->ID]);
//     Vector normalizedAvgU((centerData1->vector[0] + centerData2->vector[0]) / 2);
//     normalizedAvgU = normalizedAvgU / normalizedAvgU.length();

//     // normalized value of velocity used to not depend on initial value of velocity
//     double cos = normal * normalizedAvgU;
//     for (unsigned long nodeID : edge->innerNodeIDs) {
//         Node *node = &mesh->nodes[nodeID];
//         Data *data = &node->data;
//         if ((cos > -0.1e-9) && (cos < 0.1e-9)) {
//             assert(cell1->nodeIDToOppositeNodeID[nodeID] != -1);
//             assert(cell2->nodeIDToOppositeNodeID[nodeID] != -1);
//             double phi21 = centerData1->s1[0];
//             double phi22 = centerData2->s1[0];
//             data->s2[0] = (phi21 + phi22) / 2;
//             continue;
//         }
//         if (cos > 0) {
//             // assert(cell1->nodeIDToOppositeNodeID[nodeID] != -1);
//             Data *oppositeData1 = &mesh->nodes[cell1->nodeIDToOppositeNodeID[nodeID]].data;
//             data->s2[0] = getNewInvariantValue(
//                 data, centerData1, oppositeData1, cell1->nodeToTransferVector[nodeID], tau, true);
//         } else {
//             // assert(cell2->nodeIDToOppositeNodeID[nodeID] != -1);
//             Data *oppositeData2 = &mesh->nodes[cell2->nodeIDToOppositeNodeID[nodeID]].data;
//             data->s2[0] = getNewInvariantValue(
//                 data, centerData2, oppositeData2, cell2->nodeToTransferVector[nodeID], tau, true);
//         }
//     }
// }

// void TransferSolver::processPhase2(double tau) {
//     for (unsigned long i = 0; i < mesh->edges.size(); i++) {
//         Edge *edge = &mesh->edges[i];
//         if (edge->boundEdge) {
//             processPhase2BoundEdge(edge, tau);
//         } else {
//             processPhase2InnerEdge(edge, tau);
//         }
//     }

//     for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
//         Node *node = &mesh->nodes[i];
//         Data *data = &node->data;
//         if (!node->isApex || !node->used) {
//             continue;
//         }
//         if (node->boundNode) {
//             // TODO proper bound condition
//             data->s2[0] = 0;
//             continue;
//         }

//         double uX = 0;
//         double uY = 0;
//         for (unsigned long edgeID : node->edgeIDs) {
//             Edge *edge = &mesh->edges[edgeID];
//             long neighborInnerNodeID = edge->getNearInnerNode(node->ID);
//             Vector u = mesh->nodes[neighborInnerNodeID].data.vector[0];
//             uX += u.x;
//             uY += u.y;
//         }
//         Vector avgU(uX, uY);
//         avgU = avgU / node->edgeIDs.size();
//         Vector normalizedAvgU = avgU / avgU.length();

//         double maxCos = -1;
//         long cellIDWithMaxCos = -1;
//         for (unsigned long cellID : node->cellIDs) {
//             Vector transfer = mesh->cells[cellID].nodeToTransferVector[node->ID];
//             double cos = transfer * normalizedAvgU;
//             assert(cos >= -1 && cos <= 1);
//             if (cos >= maxCos) {
//                 maxCos = cos;
//                 cellIDWithMaxCos = cellID;
//             }
//         }
//         assert(cellIDWithMaxCos != -1);

//         Cell *cell = &mesh->cells[cellIDWithMaxCos];
//         Data *centerData = &mesh->nodes[cell->centerNodeID].data;

//         arma::vec phi0(10);
//         int pos = 0;
//         for (long nodeID : cell->nodeIDs) {
//             phi0[pos] = mesh->nodes[nodeID].data.s0[0];
//             pos++;
//         }
//         for (long edgeID : cell->edgeIDs) {
//             for (long nodeID : mesh->edges[edgeID].innerNodeIDs) {
//                 phi0[pos] = mesh->nodes[nodeID].data.s0[0];
//                 pos++;
//             }
//         }
//         phi0[pos] = centerData->s0[0];
//         arma::vec a = cell->interpolationMat * phi0;

//         double coef = -(avgU * cell->nodeToTransferVector[node->ID]) * tau;
//         Vector intersectCoords = cell->nodeToTransferVector[node->ID] * coef;
//         intersectCoords = intersectCoords + data->coords;
//         double intersectX = intersectCoords.x;
//         double intersectY = intersectCoords.y;

//         double phi2 = 0;
//         phi2 += a[0];
//         phi2 += a[1] * intersectX;
//         phi2 += a[2] * intersectY;
//         phi2 += a[3] * intersectX * intersectX;
//         phi2 += a[4] * intersectY * intersectY;
//         phi2 += a[5] * intersectX * intersectY;
//         phi2 += a[6] * intersectX * intersectX * intersectX;
//         phi2 += a[7] * intersectY * intersectY * intersectY;
//         phi2 += a[8] * intersectX * intersectY * intersectX;
//         phi2 += a[9] * intersectX * intersectY * intersectY;

//         double phiLeft0 = centerData->s0[0];
//         double phiCenter0 = centerData->s0[0];
//         double phiRight0 = data->s0[0];

//         double phiCenter1 = centerData->s1[0];

//         double h = (centerData->coords - data->coords).length();
//         double Q = (phiCenter1 - phiCenter0) / (tau / 2) +
//                    (centerData->vector[0] * cell->nodeToTransferVector[node->ID]) *
//                        (phiRight0 - phiLeft0) / h;
//         // double min = std::min(std::min(phiRight0, phiLeft0), phiCenter0) + tau * Q;
//         // double max = std::max(std::max(phiRight0, phiLeft0), phiCenter0) + tau * Q;

//         // phi2 = std::max(phi2, min);
//         // phi2 = std::min(phi2, max);

//         phi2 += tau * Q;

//         data->s2[0] = phi2;
//     }
// }

// void TransferSolver::processPhase3(double tau) {
//     for (unsigned long i = 0; i < mesh->cells.size(); i++) {
//         Cell *cell = &mesh->cells[i];
//         Node *cellCenter = &mesh->nodes[cell->centerNodeID];

//         double div = 0;
//         for (long edgeID : cell->edgeIDs) {
//             div += calcDivOnEdge(&mesh->edges[edgeID], 3) * cell->edgeToNormalDir[edgeID];
//         }
//         div /= cell->volume;

//         cellCenter->data.s2[0] = cellCenter->data.s1[0] - tau * div / 2;
//     }
// }

// void TransferSolver::prepareNextStep() {
//     for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
//         Node *node = &mesh->nodes[i];

//         Data *data = &node->data;
//         data->s0[0] = data->s2[0];
//     }
// }
