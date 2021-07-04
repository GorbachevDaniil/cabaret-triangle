#include "gtest/gtest.h"

#include "grid/mesh.hpp"

class MeshTest : public ::testing::Test {

public:
    MeshTest() : mesh_1_inner(1, false), mesh_2_inner(2, false) {
        mesh_1_inner.nodes.emplace_back(mesh_1_inner, x_1, y_1, false, true);
        mesh_1_inner.nodes.emplace_back(mesh_1_inner, x_2, y_2, false, true);
        mesh_1_inner.nodes.emplace_back(mesh_1_inner, x_3, y_3, false, true);

        mesh_1_inner.edges.emplace_back(
                mesh_1_inner, edge_id_1, mesh_1_inner.nodes[0].id, mesh_1_inner.nodes[1].id, false, mesh_1_inner.edge_inner_nodes
        );
        mesh_1_inner.edges.emplace_back(
                mesh_1_inner, edge_id_2, mesh_1_inner.nodes[1].id, mesh_1_inner.nodes[2].id, false, mesh_1_inner.edge_inner_nodes
        );
        mesh_1_inner.edges.emplace_back(
                mesh_1_inner, edge_id_3, mesh_1_inner.nodes[2].id, mesh_1_inner.nodes[0].id, false, mesh_1_inner.edge_inner_nodes
        );

        mesh_1_inner.cells.emplace_back(
                mesh_1_inner, cell_id, mesh_1_inner.nodes[0].id, mesh_1_inner.nodes[1].id, mesh_1_inner.nodes[2].id
        );

        mesh_2_inner.nodes.emplace_back(mesh_2_inner, x_1, y_1, false, true);
        mesh_2_inner.nodes.emplace_back(mesh_2_inner, x_2, y_2, false, true);
        mesh_2_inner.nodes.emplace_back(mesh_2_inner, x_3, y_3, false, true);

        mesh_2_inner.edges.emplace_back(
                mesh_2_inner, edge_id_1, mesh_2_inner.nodes[0].id, mesh_2_inner.nodes[1].id, false, mesh_2_inner.edge_inner_nodes
        );
        mesh_2_inner.edges.emplace_back(
                mesh_2_inner, edge_id_2, mesh_2_inner.nodes[1].id, mesh_2_inner.nodes[2].id, false, mesh_2_inner.edge_inner_nodes
        );
        mesh_2_inner.edges.emplace_back(
                mesh_2_inner, edge_id_3, mesh_2_inner.nodes[2].id, mesh_2_inner.nodes[0].id, false, mesh_2_inner.edge_inner_nodes
        );

        mesh_2_inner.cells.emplace_back(
                mesh_2_inner, cell_id, mesh_2_inner.nodes[0].id, mesh_2_inner.nodes[1].id, mesh_2_inner.nodes[2].id
        );
    }

protected:
    Mesh mesh_1_inner;
    Mesh mesh_2_inner;
    const double x_1 = 0;
    const double y_1 = 0;
    const double x_2 = 1;
    const double y_2 = 0;
    const double x_3 = 0;
    const double y_3 = 1;

    const long edge_id_1 = 0;
    const long edge_id_2 = 1;
    const long edge_id_3 = 2;
    const long cell_id = 0;
};

TEST_F(MeshTest, return_next_id_in_a_circular_manner) {
    EXPECT_EQ(mesh_1_inner.cells[0].get_next_node_id(0), mesh_1_inner.nodes[1].id);
    EXPECT_EQ(mesh_1_inner.cells[0].get_next_node_id(1), mesh_1_inner.nodes[2].id);
    EXPECT_EQ(mesh_1_inner.cells[0].get_next_node_id(2), mesh_1_inner.nodes[0].id);
}

TEST_F(MeshTest, return_prev_id_in_a_circular_manner) {
    EXPECT_EQ(mesh_1_inner.cells[0].get_prev_node_id(0), mesh_1_inner.nodes[2].id);
    EXPECT_EQ(mesh_1_inner.cells[0].get_prev_node_id(1), mesh_1_inner.nodes[0].id);
    EXPECT_EQ(mesh_1_inner.cells[0].get_prev_node_id(2), mesh_1_inner.nodes[1].id);
}

TEST_F(MeshTest, die_if_next_id_more_than_size) {
    ASSERT_DEATH(mesh_1_inner.cells[0].get_next_node_id(mesh_1_inner.cells[0].node_ids.size()), "");
}

TEST_F(MeshTest, die_if_next_id_negative) {
    ASSERT_DEATH(mesh_1_inner.cells[0].get_next_node_id(-1), "");
}

TEST_F(MeshTest, die_if_prev_id_more_than_size) {
    ASSERT_DEATH(mesh_1_inner.cells[0].get_prev_node_id(mesh_1_inner.cells[0].node_ids.size()), "");
}

TEST_F(MeshTest, die_if_prev_id_negative) {
    ASSERT_DEATH(mesh_1_inner.cells[0].get_prev_node_id(-1), "");
}

TEST_F(MeshTest, get_edge_ordered_node_ids_positive_case) {
    Cell& cell = mesh_1_inner.cells[0];
    Edge& edge_1 = mesh_1_inner.edges[edge_id_1];
    Edge& edge_2 = mesh_1_inner.edges[edge_id_2];
    Edge& edge_3 = mesh_1_inner.edges[edge_id_3];

    std::vector<long> forward_node_ids_edge_1{{edge_1.end_node_ids[0], edge_1.end_node_ids[1]}};
    std::vector<long> backward_node_ids_edge_1{{edge_1.end_node_ids[1], edge_1.end_node_ids[0]}};

    std::vector<long> ordered_node_ids_edge_1_from_forward = cell.get_edge_ordered_node_ids(forward_node_ids_edge_1);
    EXPECT_EQ(ordered_node_ids_edge_1_from_forward[0], mesh_1_inner.nodes[0].id);
    EXPECT_EQ(ordered_node_ids_edge_1_from_forward[1], mesh_1_inner.nodes[1].id);

    std::vector<long> ordered_node_ids_edge_1_from_backward = cell.get_edge_ordered_node_ids(backward_node_ids_edge_1);
    EXPECT_EQ(ordered_node_ids_edge_1_from_backward[0], mesh_1_inner.nodes[0].id);
    EXPECT_EQ(ordered_node_ids_edge_1_from_backward[1], mesh_1_inner.nodes[1].id);

    std::vector<long> forward_node_ids_edge_2{{edge_2.end_node_ids[0], edge_2.end_node_ids[1]}};
    std::vector<long> backward_node_ids_edge_2{{edge_2.end_node_ids[1], edge_2.end_node_ids[0]}};

    std::vector<long> ordered_node_ids_edge_2_from_forward = cell.get_edge_ordered_node_ids(forward_node_ids_edge_2);
    EXPECT_EQ(ordered_node_ids_edge_2_from_forward[0], mesh_1_inner.nodes[1].id);
    EXPECT_EQ(ordered_node_ids_edge_2_from_forward[1], mesh_1_inner.nodes[2].id);

    std::vector<long> ordered_node_ids_edge_2_from_backward = cell.get_edge_ordered_node_ids(backward_node_ids_edge_2);
    EXPECT_EQ(ordered_node_ids_edge_2_from_backward[0], mesh_1_inner.nodes[1].id);
    EXPECT_EQ(ordered_node_ids_edge_2_from_backward[1], mesh_1_inner.nodes[2].id);

    std::vector<long> forward_node_ids_edge_3{{edge_3.end_node_ids[0], edge_3.end_node_ids[1]}};
    std::vector<long> backward_node_ids_edge_3{{edge_3.end_node_ids[1], edge_3.end_node_ids[0]}};

    std::vector<long> ordered_node_ids_edge_3_from_forward = cell.get_edge_ordered_node_ids(forward_node_ids_edge_3);
    EXPECT_EQ(ordered_node_ids_edge_3_from_forward[0], mesh_1_inner.nodes[2].id);
    EXPECT_EQ(ordered_node_ids_edge_3_from_forward[1], mesh_1_inner.nodes[0].id);

    std::vector<long> ordered_node_ids_edge_3_from_backward = cell.get_edge_ordered_node_ids(backward_node_ids_edge_3);
    EXPECT_EQ(ordered_node_ids_edge_3_from_backward[0], mesh_1_inner.nodes[2].id);
    EXPECT_EQ(ordered_node_ids_edge_3_from_backward[1], mesh_1_inner.nodes[0].id);
}

TEST_F(MeshTest, get_edge_ordered_node_ids_negative_case) {
    Cell& cell = mesh_1_inner.cells[0];

    ASSERT_DEATH(cell.get_edge_ordered_node_ids({0}), "");
    ASSERT_DEATH(cell.get_edge_ordered_node_ids({0, 1, 2}), "");
    ASSERT_DEATH(cell.get_edge_ordered_node_ids({0, 0}), "");
    ASSERT_DEATH(cell.get_edge_ordered_node_ids({0, -1}), "");
}

TEST_F(MeshTest, one_inner_node_doesnt_have_opposite) {
    Cell& cell = mesh_1_inner.cells[0];
    Edge& edge_1 = mesh_1_inner.edges[edge_id_1];
    Edge& edge_2 = mesh_1_inner.edges[edge_id_2];
    Edge& edge_3 = mesh_1_inner.edges[edge_id_3];

    EXPECT_EQ(cell.node_id_to_opposite_node_id[edge_1.inner_node_ids[0]], -1);
    EXPECT_EQ(cell.node_id_to_opposite_node_id[edge_2.inner_node_ids[0]], -1);
    EXPECT_EQ(cell.node_id_to_opposite_node_id[edge_3.inner_node_ids[0]], -1);
}

TEST_F(MeshTest, two_inner_node_have_opposite) {
    Cell& cell = mesh_2_inner.cells[0];
    Edge& edge_1 = mesh_2_inner.edges[edge_id_1];
    Edge& edge_2 = mesh_2_inner.edges[edge_id_2];
    Edge& edge_3 = mesh_2_inner.edges[edge_id_3];

    EXPECT_EQ(6, cell.node_id_to_opposite_node_id[edge_1.inner_node_ids[0]]);
    EXPECT_EQ(7, cell.node_id_to_opposite_node_id[edge_1.inner_node_ids[1]]);
    EXPECT_EQ(8, cell.node_id_to_opposite_node_id[edge_2.inner_node_ids[0]]);
    EXPECT_EQ(3, cell.node_id_to_opposite_node_id[edge_2.inner_node_ids[1]]);
    EXPECT_EQ(4, cell.node_id_to_opposite_node_id[edge_3.inner_node_ids[0]]);
    EXPECT_EQ(5, cell.node_id_to_opposite_node_id[edge_3.inner_node_ids[1]]);
}

TEST_F(MeshTest, edge_with_one_inner_node_must_have_three_nodes_in_total) {
    Edge& edge_1 = mesh_1_inner.edges[edge_id_1];
    Edge& edge_2 = mesh_1_inner.edges[edge_id_2];
    Edge& edge_3 = mesh_1_inner.edges[edge_id_3];

    EXPECT_EQ(3, edge_1.node_ids.size());
    EXPECT_EQ(3, edge_2.node_ids.size());
    EXPECT_EQ(3, edge_3.node_ids.size());

    EXPECT_EQ(1, edge_1.inner_node_ids.size());
    EXPECT_EQ(1, edge_2.inner_node_ids.size());
    EXPECT_EQ(1, edge_3.inner_node_ids.size());

    EXPECT_EQ(2, edge_1.end_node_ids.size());
    EXPECT_EQ(2, edge_2.end_node_ids.size());
    EXPECT_EQ(2, edge_3.end_node_ids.size());
}

TEST_F(MeshTest, edge_with_two_inner_node_must_have_four_nodes_in_total) {
    Edge& edge_1 = mesh_2_inner.edges[edge_id_1];
    Edge& edge_2 = mesh_2_inner.edges[edge_id_2];
    Edge& edge_3 = mesh_2_inner.edges[edge_id_3];

    EXPECT_EQ(4, edge_1.node_ids.size());
    EXPECT_EQ(4, edge_2.node_ids.size());
    EXPECT_EQ(4, edge_3.node_ids.size());

    EXPECT_EQ(2, edge_1.inner_node_ids.size());
    EXPECT_EQ(2, edge_2.inner_node_ids.size());
    EXPECT_EQ(2, edge_3.inner_node_ids.size());

    EXPECT_EQ(2, edge_1.end_node_ids.size());
    EXPECT_EQ(2, edge_2.end_node_ids.size());
    EXPECT_EQ(2, edge_3.end_node_ids.size());
}

TEST_F(MeshTest, calculate_edge_normals) {
    mesh_1_inner.calculate_edges_normals();

    EXPECT_EQ(mesh_1_inner.cells[0].edge_to_normal[edge_id_1].x, 0);
    EXPECT_EQ(mesh_1_inner.cells[cell_id].edge_to_normal[edge_id_1].y, -1);
    EXPECT_EQ(mesh_1_inner.cells[cell_id].edge_to_normal[edge_id_2].x, 1 / sqrt(2));
    EXPECT_EQ(mesh_1_inner.cells[cell_id].edge_to_normal[edge_id_2].y, 1 / sqrt(2));
    EXPECT_EQ(mesh_1_inner.cells[cell_id].edge_to_normal[edge_id_3].x, -1);
    EXPECT_EQ(mesh_1_inner.cells[cell_id].edge_to_normal[edge_id_3].y, 0);
}

TEST_F(MeshTest, calculate_transfer_vectors) {
    mesh_1_inner.calculate_transfer_vectors();

    double x_center = (x_1 + x_2 + x_3) / 3;
    double y_center = (y_1 + y_2 + y_3) / 3;
    double transfer_x = x_center - (x_1 + x_2) / 2;
    double transfer_y = y_center - (y_1 + y_2) / 2;
    double transfer_len = Vector::length(transfer_x, transfer_y);
    EXPECT_EQ(mesh_1_inner.cells[cell_id].node_to_transfer_vector[mesh_1_inner.nodes[3].id].x, transfer_x / transfer_len);
    EXPECT_EQ(mesh_1_inner.cells[cell_id].node_to_transfer_vector[mesh_1_inner.nodes[3].id].y, transfer_y / transfer_len);

    transfer_x = x_center - (x_2 + x_3) / 2;
    transfer_y = y_center - (y_2 + y_3) / 2;
    transfer_len = Vector::length(transfer_x, transfer_y);
    EXPECT_EQ(mesh_1_inner.cells[cell_id].node_to_transfer_vector[mesh_1_inner.nodes[4].id].x, transfer_x / transfer_len);
    EXPECT_EQ(mesh_1_inner.cells[cell_id].node_to_transfer_vector[mesh_1_inner.nodes[4].id].y, transfer_y / transfer_len);

    transfer_x = x_center - (x_1 + x_3) / 2;
    transfer_y = y_center - (y_1 + y_3) / 2;
    transfer_len = Vector::length(transfer_x, transfer_y);
    EXPECT_EQ(mesh_1_inner.cells[cell_id].node_to_transfer_vector[mesh_1_inner.nodes[5].id].x, transfer_x / transfer_len);
    EXPECT_EQ(mesh_1_inner.cells[cell_id].node_to_transfer_vector[mesh_1_inner.nodes[5].id].y, transfer_y / transfer_len);
}
