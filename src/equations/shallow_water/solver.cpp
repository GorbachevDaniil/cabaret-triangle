#include "equations/shallow_water/solver.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

void ShallowWaterSolver::calc_tau() {
    double tau = std::numeric_limits<double>::max();
    for (auto& cell : mesh_.cells) {
        long center_node_id = cell.center_node_id;
        Node* center_node = &mesh_.nodes[center_node_id];
        for (long edge_id : cell.edge_ids) {
            for (long inner_node_id : mesh_.edges[edge_id].inner_node_ids) {
                Node* node = &mesh_.nodes[inner_node_id];

                double h = 2 * Vector::length(node->coords.x - center_node->coords.x,
                                              node->coords.y - center_node->coords.y);

                double lambda_r = calc_lambda_r(mesh_.s0[center_node_id][0],
                                                mesh_.v0[center_node_id][0],
                                                cell.node_to_transfer_vector[inner_node_id]);
                double node_tau_r = h / std::abs(lambda_r);
                if (tau > node_tau_r) {
                    tau = node_tau_r;
                }

                double lambda_q = calc_lambda_q(mesh_.s0[center_node_id][0],
                                                mesh_.v0[center_node_id][0],
                                                cell.node_to_transfer_vector[inner_node_id]);
                double node_tau_q = h / std::abs(lambda_q);
                if (tau > node_tau_q) {
                    tau = node_tau_q;
                }
            }
        }
    }
    assert(tau != std::numeric_limits<double>::max());
    tau_ = cfl_ * tau;
}

double ShallowWaterSolver::calc_div_1(double h, double u_x, double u_y, Vector div_coef) const {
    return h * (u_x * div_coef.x + u_y * div_coef.y);
}

double ShallowWaterSolver::calc_div_2(double h, double u_x, double u_y, Vector div_coef) const {
    double f_x = pow(u_x, 2) + g_ * h / 2;
    double f_y = u_x * u_y;
    return h * (f_x * div_coef.x + f_y * div_coef.y);
}

double ShallowWaterSolver::calc_div_3(double h, double u_x, double u_y, Vector div_coef) const {
    double f_x = u_x * u_y;
    double f_y = pow(u_y, 2) + g_ * h / 2;
    return h * (f_x * div_coef.x + f_y * div_coef.y);
}

void ShallowWaterSolver::process_phase_1() {
    for (auto& cell : mesh_.cells) {
        double div1 = 0;
        double div2 = 0;
        double div3 = 0;
        for (long edge_id : cell.edge_ids) {
            Edge *edge = &mesh_.edges[edge_id];

            double h = 0;
            double u_x = 0;
            double u_y = 0;

            for (unsigned long used_node_id : edge->used_node_ids) {
                h += mesh_.s0[used_node_id][0];
                u_x += mesh_.v0[used_node_id][0].x;
                u_y += mesh_.v0[used_node_id][0].y;
            }
            h /= (double) edge->used_node_ids.size();
            u_x /= (double) edge->used_node_ids.size();
            u_y /= (double) edge->used_node_ids.size();

            Vector div_coef = cell.edge_to_div_coef[edge_id];
            div1 += calc_div_1(h, u_x, u_y, div_coef);
            div2 += calc_div_2(h, u_x, u_y, div_coef);
            div3 += calc_div_3(h, u_x, u_y, div_coef);
        }
        div1 /= cell.volume;
        div2 /= cell.volume;
        div3 /= cell.volume;

        long center_node_id = cell.center_node_id;

        double h = mesh_.s0[center_node_id][0];
        double u_x = mesh_.v0[center_node_id][0].x;
        double u_y = mesh_.v0[center_node_id][0].y;

        double new_h = h - tau_ * div1 / 2;
        double new_u_x = (h * u_x - tau_ * div2 / 2) / new_h;
        double new_u_y = (h * u_y - tau_ * div3 / 2) / new_h;

        mesh_.s1[center_node_id][0] = new_h;
        mesh_.v1[center_node_id][0].x = new_u_x;
        mesh_.v1[center_node_id][0].y = new_u_y;
    }
}

double ShallowWaterSolver::calc_inv_r(double G, double h, const Vector& u, const Vector& n) const {
    return (u.x * n.x + u.y * n.y) + G * h;
}

double ShallowWaterSolver::calc_inv_q(double G, double h, const Vector& u, const Vector& n) const {
    return (u.x * n.x + u.y * n.y) - G * h;
}

double ShallowWaterSolver::calc_inv_s(const Vector& u, const Vector& n) const {
    return -u.x * n.y + u.y * n.x;
}

double ShallowWaterSolver::calc_lambda_r(double h, const Vector& u, const Vector& n) const {
    return (u.x * n.x + u.y * n.y) + sqrt(g_ * h);
}

double ShallowWaterSolver::calc_lambda_q(double h, const Vector& u, const Vector& n) const {
    return (u.x * n.x + u.y * n.y) - sqrt(g_ * h);
}

double ShallowWaterSolver::calc_lambda_s(const Vector& u, const Vector& n) const {
    return (u.x * n.x + u.y * n.y);
}

void ShallowWaterSolver::process_phase_2_bound(const Edge& edge) {
    for (const long node_id : edge.inner_node_ids) {
        mesh_.v2[node_id][0].x = 0;
        mesh_.v2[node_id][0].y = 0;
        mesh_.s2[node_id][0] = 1;
    }
}

std::array<std::array<double, 3>, 2> ShallowWaterSolver::calc_invs(const double& G, const Vector& n,
                                                                   Cell& cell, long node_id) {
    double pan_coef = 0.1;
    std::array<std::array<double, 3>, 2> array = {{
        {0, 0, 0},
        {0, 0, 0}
    }};

    long center_node_id = cell.center_node_id;
    long opposite_node_id = cell.node_id_to_opposite_node_id[node_id];

    double lambda_r = calc_lambda_r(mesh_.s1[center_node_id][0], mesh_.v1[center_node_id][0], n);
    double inv_r_opposite = calc_inv_r(G, mesh_.s0[opposite_node_id][0], mesh_.v0[opposite_node_id][0], n);
    double inv_r_center = calc_inv_r(G, mesh_.s1[center_node_id][0], mesh_.v1[center_node_id][0], n);
    double inv_r_node = (2 * inv_r_center - (1 - pan_coef) * inv_r_opposite) / (1 + pan_coef);
    array[0][0] = lambda_r;
    array[1][0] = inv_r_node;

    double lambda_q = calc_lambda_q(mesh_.s1[center_node_id][0], mesh_.v1[center_node_id][0], n);
    double inv_q_opposite = calc_inv_q(G, mesh_.s0[opposite_node_id][0], mesh_.v0[opposite_node_id][0], n);
    double inv_q_center = calc_inv_q(G, mesh_.s1[center_node_id][0], mesh_.v1[center_node_id][0], n);
    double inv_q_node = (2 * inv_q_center - (1 - pan_coef) * inv_q_opposite) / (1 + pan_coef);
    array[0][1] = lambda_q;
    array[1][1] = inv_q_node;

    double lambda_s = calc_lambda_s(mesh_.v1[center_node_id][0], n);
    double inv_s_opposite = calc_inv_s(mesh_.v0[opposite_node_id][0], n);
    double inv_s_center = calc_inv_s(mesh_.v1[center_node_id][0], n);
    double inv_s_node = (2 * inv_s_center - (1 - pan_coef) * inv_s_opposite) / (1 + pan_coef);
    array[0][2] = lambda_s;
    array[1][2] = inv_s_node;

    return array;
}

void ShallowWaterSolver::process_phase_2_inner(const Edge& edge) {
    Vector n;
    Cell& cell_1 = mesh_.cells[edge.cell_ids[0]];
    Cell& cell_2 = mesh_.cells[edge.cell_ids[1]];
    n.x = (cell_1.edge_to_normal[edge.id].x - cell_2.edge_to_normal[edge.id].x) / 2;
    n.y = (cell_1.edge_to_normal[edge.id].y - cell_2.edge_to_normal[edge.id].y) / 2;

    double G_1 = sqrt(g_ / mesh_.s1[cell_1.center_node_id][0]);
    double G_2 = sqrt(g_ / mesh_.s1[cell_2.center_node_id][0]);

    for (const long node_id : edge.inner_node_ids) {
        std::array<std::array<double, 3>, 2> invs_1 = calc_invs(G_1, n, cell_1, node_id);
        std::array<std::array<double, 3>, 2> invs_2 = calc_invs(G_2, n, cell_2, node_id);

        double G_L;
        double G_R;

        double R;
        if (invs_1[0][0] > 0) {
            R = invs_1[1][0];
            G_L = G_1;
        } else {
            R = invs_2[1][0];
            G_L = G_2;
        }

        double Q;
        if (invs_1[0][1] > 0) {
            Q = invs_1[1][1];
            G_R = G_1;
        } else {
            Q = invs_2[1][1];
            G_R = G_2;
        }

        double S;
        if (invs_1[0][2] > 0) {
            S = invs_1[1][2];
        } else {
            S = invs_2[1][2];
        }

        double coef = (G_R * R + G_L * Q) / (G_R + G_L);
        mesh_.v2[node_id][0].x = n.x * coef - n.y * S;
        mesh_.v2[node_id][0].y = n.y * coef + n.x * S;
        mesh_.s2[node_id][0] = (R-Q) / (G_R + G_L);
    }

//    Cell& cell_1 = mesh_.cells[edge.cell_ids[0]];
//    Cell& cell_2 = mesh_.cells[edge.cell_ids[1]];
//
//    double G_1 = sqrt(g_ / mesh_.s1[cell_1.center_node_id][0]);
//    double G_2 = sqrt(g_ / mesh_.s1[cell_2.center_node_id][0]);
//
//    for (const long node_id : edge.inner_node_ids) {
//        Vector n_1;
//        n_1.x = cell_1.node_to_transfer_vector[node_id].x;
//        n_1.y = cell_1.node_to_transfer_vector[node_id].y;
//        std::array<std::array<double, 3>, 2> invs_1 = calc_invs(G_1, n_1, cell_1, node_id);
//
//        Vector n_2;
//        n_2.x = cell_2.node_to_transfer_vector[node_id].x;
//        n_2.y = cell_2.node_to_transfer_vector[node_id].y;
//        std::array<std::array<double, 3>, 2> invs_2 = calc_invs(G_2, n_2, cell_2, node_id);
//
//        if ((invs_1[0][0] == invs_2[0][0]) &&
//            (invs_1[0][1] == invs_2[0][1]) &&
//            (invs_1[0][2] == invs_2[0][2])) {
//            mesh_.v2[node_id][0].x = (mesh_.v1[cell_1.center_node_id][0].x + mesh_.v1[cell_1.center_node_id][0].x) / 2;
//            mesh_.v2[node_id][0].y = (mesh_.v1[cell_1.center_node_id][0].y + mesh_.v1[cell_1.center_node_id][0].y) / 2;
//            mesh_.s2[node_id][0] = (mesh_.s1[cell_1.center_node_id][0] + mesh_.s1[cell_1.center_node_id][0]) / 2;
//            continue;
//        }
//
//        double G_R;
//        double G_Q;
//
//        double R;
//        Vector n_R;
//        if (invs_1[0][0] > 0) {
//            R = invs_1[1][0];
//            G_R = G_1;
//            n_R.x = n_1.x;
//            n_R.y = n_1.y;
//        } else {
//            R = invs_2[1][0];
//            G_R = G_2;
//            n_R.x = n_2.x;
//            n_R.y = n_2.y;
//        }
//
//        double Q;
//        Vector n_Q;
//        if (invs_1[0][1] > 0) {
//            Q = invs_1[1][1];
//            G_Q = G_1;
//            n_Q.x = n_1.x;
//            n_Q.y = n_1.y;
//        } else {
//            Q = invs_2[1][1];
//            G_Q = G_2;
//            n_Q.x = n_2.x;
//            n_Q.y = n_2.y;
//        }
//
//
//        double S;
//        Vector n_S;
//        if (invs_1[0][2] > 0) {
//            S = invs_1[1][2];
//            n_S.x = n_1.x;
//            n_S.y = n_1.y;
//        } else {
//            S = invs_2[1][2];
//            n_S.x = n_2.x;
//            n_S.y = n_2.y;
//        }
//
//        double coef = G_R * (n_Q.x * n_S.x + n_Q.y * n_S.y) +
//                      G_Q * (n_R.x * n_S.x + n_R.y * n_S.y);
//        mesh_.v2[node_id][0].x = (G_R * (n_S.x * Q - n_Q.y * S) +
//                                  G_Q * (n_S.x * R - n_R.y * S)) / coef;
//        mesh_.v2[node_id][0].y = (G_R * (n_S.y * Q + n_Q.x * S) +
//                                  G_Q * (n_S.y * R + n_R.x * S)) / coef;
//        mesh_.s2[node_id][0] = (n_S.x * (n_Q.x * R - n_R.x * Q) +
//                                n_S.y * (n_Q.y * R - n_R.y * Q) +
//                                S * (n_R.x * n_Q.y - n_R.y * n_Q.x)) / coef;
//    }
}

void ShallowWaterSolver::process_phase_2() {
    for (unsigned long i = 0; i < mesh_.edges.size(); i++) {
        Edge& edge = mesh_.edges[i];
        if (!edge.is_bound) {
            process_phase_2_inner(edge);
        } else {
            process_phase_2_bound(edge);
        }
    }
}

void ShallowWaterSolver::process_phase_3() {
    for (auto& cell : mesh_.cells) {
        double div1 = 0;
        double div2 = 0;
        double div3 = 0;
        for (long edge_id : cell.edge_ids) {
            Edge *edge = &mesh_.edges[edge_id];

            double h = 0;
            double u_x = 0;
            double u_y = 0;

            for (unsigned long used_node_id : edge->used_node_ids) {
                h += mesh_.s1[used_node_id][0];
                u_x += mesh_.v1[used_node_id][0].x;
                u_y += mesh_.v1[used_node_id][0].y;
            }
            h /= (double) edge->used_node_ids.size();
            u_x /= (double) edge->used_node_ids.size();
            u_y /= (double) edge->used_node_ids.size();

            Vector div_coef = cell.edge_to_div_coef[edge_id];
            div1 += calc_div_1(h, u_x, u_y, div_coef);
            div2 += calc_div_2(h, u_x, u_y, div_coef);
            div3 += calc_div_3(h, u_x, u_y, div_coef);
        }
        div1 /= cell.volume;
        div2 /= cell.volume;
        div3 /= cell.volume;

        long center_node_id = cell.center_node_id;

        double h = mesh_.s1[center_node_id][0];
        double u_x = mesh_.v1[center_node_id][0].x;
        double u_y = mesh_.v1[center_node_id][0].y;

        double new_h = h - tau_ * div1 / 2;
        double new_u_x = (h * u_x - tau_ * div2 / 2) / new_h;
        double new_u_y = (h * u_y - tau_ * div3 / 2) / new_h;

        mesh_.s2[center_node_id][0] = new_h;
        mesh_.v2[center_node_id][0].x = new_u_x;
        mesh_.v2[center_node_id][0].y = new_u_y;
    }
}

void ShallowWaterSolver::prepare_next_step() {
    for (unsigned long i = 0; i < mesh_.nodes.size(); i++) {
        for (unsigned long j = 0; j < mesh_.s0[i].size(); j++) {
            mesh_.s0[i][j] = mesh_.s2[i][j];
            mesh_.s1[i][j] = 0;
            mesh_.s2[i][j] = 0;
        }
        for (unsigned long j = 0; j < mesh_.v0[i].size(); j++) {
            mesh_.v0[i][j].x = mesh_.v2[i][j].x;
            mesh_.v0[i][j].y = mesh_.v2[i][j].y;
            mesh_.v1[i][j].x = 0;
            mesh_.v1[i][j].y = 0;
            mesh_.v2[i][j].x = 0;
            mesh_.v2[i][j].y = 0;
        }
    }
}