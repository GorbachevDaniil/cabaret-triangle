#include <chrono>
#include <iostream>

#include "grid/mesh.hpp"
#include "equations/shallow_water/initializer.hpp"
#include "equations/shallow_water/output.hpp"
#include "equations/shallow_water/solver.hpp"

int main() {
    std::cout.precision(30);

    int steps = 500;
    double time_to_stop = 0.0;

    int write_period = 1;
    bool write_conservative = true;
    bool write_flux = true;

    int edge_inner_nodes = 2;
    bool apex_nodes_used = false;

    Mesh mesh{edge_inner_nodes, apex_nodes_used};
    mesh.init_mesh(
        "/Users/gorbdan/workspace/projects/cabaret/cabaret-triangle/resources/mesh/regular/square -1x1/4080/Mesh.node",
        "/Users/gorbdan/workspace/projects/cabaret/cabaret-triangle/resources/mesh/regular/square -1x1/4080/Mesh.edge",
        "/Users/gorbdan/workspace/projects/cabaret/cabaret-triangle/resources/mesh/regular/square -1x1/4080/Mesh.ele"
    );
    std::cout << "number of cells = " << mesh.cells.size() << std::endl;

    mesh.calculate_edges_normals();
    mesh.calculate_transfer_vectors();

    int task = 2;
    switch (task) {
        case 2: {
            ShallowWaterInitializer initializer = ShallowWaterInitializer();
            initializer.initialize(mesh);

            ShallowWaterOutput output = ShallowWaterOutput{write_conservative, write_flux};
            output.write_paraview(mesh, 0, 0);

            double cfl = 0.3;
            double g = 1.0;
            ShallowWaterSolver solver = ShallowWaterSolver(mesh, cfl, g);

            std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

            std::chrono::duration<double> durCalcTau = std::chrono::duration<double>(0);
            std::chrono::duration<double> durPhase1 = std::chrono::duration<double>(0);
            std::chrono::duration<double> durPhase2 = std::chrono::duration<double>(0);
            std::chrono::duration<double> durPhase3 = std::chrono::duration<double>(0);

            double tau;
            double time = 0;
            int real_steps = 0;
            for (int i = 0; i < steps; i++) {
                std::chrono::system_clock::time_point startCalcTau = std::chrono::system_clock::now();
                solver.calc_tau();
                tau = solver.get_tau();
                std::cout << "step = " << i + 1 << " tau = " << tau << " time = " << time << std::endl;
                time += tau;
                std::chrono::system_clock::time_point stopCalcTau = std::chrono::system_clock::now();
                durCalcTau +=
                    std::chrono::duration_cast<std::chrono::duration<double>>(stopCalcTau - startCalcTau);

                std::chrono::system_clock::time_point startPhase1 = std::chrono::system_clock::now();
                solver.process_phase_1();
                std::chrono::system_clock::time_point stopPhase1 = std::chrono::system_clock::now();
                durPhase1 +=
                    std::chrono::duration_cast<std::chrono::duration<double>>(stopPhase1 - startPhase1);

                std::chrono::system_clock::time_point startPhase2 = std::chrono::system_clock::now();
                solver.process_phase_2();
                std::chrono::system_clock::time_point stopPhase2 = std::chrono::system_clock::now();
                durPhase2 +=
                    std::chrono::duration_cast<std::chrono::duration<double>>(stopPhase2 - startPhase2);

                std::chrono::system_clock::time_point startPhase3 = std::chrono::system_clock::now();
                solver.process_phase_3();
                std::chrono::system_clock::time_point stopPhase3 = std::chrono::system_clock::now();
                durPhase3 +=
                    std::chrono::duration_cast<std::chrono::duration<double>>(stopPhase3 - startPhase3);

                solver.prepare_next_step();

                real_steps += 1;

                if ((i + 1) % write_period == 0) {
                    output.write_paraview(mesh, time, i + 1);
                }
                if (time_to_stop > 0 && time >= time_to_stop) {
                    break;
                }
            }
            std::chrono::system_clock::time_point stop = std::chrono::system_clock::now();
            std::chrono::duration<double> dur = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);

            std::cout << "it took " << dur.count() / real_steps << " sec to process one step" << std::endl;
            std::cout << "it took " << durCalcTau.count() / real_steps << " sec to calculate tau" << std::endl;
            std::cout << "it took " << durPhase1.count() / real_steps << " sec to process phase 1" << std::endl;
            std::cout << "it took " << durPhase2.count() / real_steps << " sec to process phase 2" << std::endl;
            std::cout << "it took " << durPhase3.count() / real_steps << " sec to process phase 3" << std::endl;
            break;
        }
        default: {
            std::cout << "such task not supported" << std::endl;
            break;
        }
    }

    return 0;
}
